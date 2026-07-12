"""Dictate a question to the copilot - speech-to-text with faster-whisper, on the server's CPU.

App-layer only: `analysis/` stays UI-agnostic and never imports fastapi, so the transcription for
the copilot's ask box lives here in `backend/`, next to `app.py`.

**Optional dependency.** `available()` is False when `faster_whisper` isn't importable; `POST /api/stt`
then 503s and the UI hides the mic button. Install it into the env that serves `/api`:

    uv pip install faster-whisper          # 7 packages, NO torch (it runs on CTranslate2)

Model weights are pulled from HuggingFace on first use into `$HF_HOME` - point that at a volume with
room, never a quota'd home directory.

Audio never leaves the server: the browser posts the recording, we transcribe locally and return
text. That is the whole reason this is not the browser's built-in `SpeechRecognition`, which in
Chrome uploads the audio to Google - unacceptable for spoken questions about patient sections.

Measured on a compute node (int8, 4 threads), on a 2.5 s question:
  tiny.en   0.36 s   "Are the T-cells excluded from the tumor core?"   (hyphenates)
  base.en   0.41 s   "Are the T cells excluded from the tumor core?"   <- default
  small.en  1.22 s   identical to base.en, 3x the cost
`base.en` it is. Override with SPATIALSCRIBE_STT_MODEL (e.g. `small`, `medium.en`, or a multilingual
`base` if you want anything other than English).

Env:
    SPATIALSCRIBE_STT_MODEL    faster-whisper model id (default "base.en")
    SPATIALSCRIBE_STT_THREADS  CTranslate2 CPU threads (default 4)
"""
from __future__ import annotations

import io
import os
import threading
from functools import lru_cache

MAX_BYTES = 12 * 1024 * 1024   # ~10 min of opus; a spoken question is a few tens of KB
_LOCK = threading.Lock()       # one model, shared - serialize inference
_MODEL = None


class RecordingTooLarge(ValueError):
    """Raised only for the size guard.

    Needed because PyAV raises ValueError *subclasses* (`av.error.InvalidDataError`) on an
    undecodable blob: a bare `except ValueError` in the endpoint reported 9 bytes of garbage as
    "recording too large" (HTTP 413), sending the user to chase a size limit that was never the
    problem. Catch this first, let every other failure fall through as a decode error.
    """


@lru_cache(maxsize=1)
def available() -> bool:
    """True when faster-whisper is installed. Cheap - no model import, no weights download."""
    import importlib.util
    return importlib.util.find_spec("faster_whisper") is not None


def _model():
    """Lazily load (and keep) the model. First call downloads weights and costs ~5 s."""
    global _MODEL
    if _MODEL is None:
        from faster_whisper import WhisperModel
        _MODEL = WhisperModel(
            os.environ.get("SPATIALSCRIBE_STT_MODEL", "base.en"),
            device="cpu", compute_type="int8",
            cpu_threads=int(os.environ.get("SPATIALSCRIBE_STT_THREADS", "4")),
        )
    return _MODEL


def warm() -> None:
    """Load the model on a background daemon thread, so the first dictated question doesn't wait
    ~5 s for weights. Best-effort: a missing package or a failed download never blocks startup."""
    if not available():
        return

    def _go():
        try:
            with _LOCK:
                _model()
        except Exception:
            pass
    threading.Thread(target=_go, daemon=True).start()


def transcribe(audio: bytes) -> str:
    """Transcribe a browser recording to text.

    ``audio`` is whatever ``MediaRecorder`` produced - in practice webm/opus, sometimes ogg/opus or
    mp4/aac depending on the browser. faster-whisper decodes it in-process via PyAV, so there is **no
    ffmpeg binary** to install (there is none on the nodes) and no temp file to clean up.

    Returns the stripped transcript, or "" for silence. Never raises on an empty result - a user who
    taps the mic and says nothing should get an empty box back, not a 500.
    """
    if not audio:
        return ""
    if len(audio) > MAX_BYTES:
        raise RecordingTooLarge(f"recording too large: {len(audio)} bytes (max {MAX_BYTES})")
    with _LOCK:
        # beam_size=1 (greedy): on a short, clean, close-mic question it matches beam search and is
        # ~2x faster. vad_filter drops the leading/trailing silence around the click-to-talk.
        segments, _info = _model().transcribe(io.BytesIO(audio), beam_size=1, vad_filter=True)
        return " ".join(s.text.strip() for s in segments).strip()


if __name__ == "__main__":   # self-check: python backend/stt.py <audio-file>
    import sys
    assert available(), "faster-whisper not installed - see the module docstring"
    if len(sys.argv) < 2:
        raise SystemExit("usage: python backend/stt.py <wav|webm|ogg|mp3>")
    text = transcribe(open(sys.argv[1], "rb").read())
    assert text, "transcribed to nothing - is the file silent?"
    print(f"ok - {text!r}")
