"""Speech-to-text for the copilot's ask box (backend/stt.py).

The model itself is an optional dep and a ~5 s load, so the transcription tests are skipped unless
faster-whisper is installed. What is checked without it: the availability gate, the size guard, and
that silence/empty input degrade to "" rather than raising (a user who taps the mic and says nothing
must get an empty box, not a 500).

The transcription tests synthesize their own audio, so there is no binary fixture in the repo. The
webm/opus one is the load-bearing case: that is what a browser MediaRecorder actually posts, and it
must decode in-process via PyAV, since no the cluster node has an ffmpeg binary.
"""
from __future__ import annotations

import io
import math
import struct
import sys
import wave
from pathlib import Path

import pytest

# backend/ lives at the repo root (a sibling of tests/), not under src - put it on the path.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from backend import stt  # noqa: E402

needs_model = pytest.mark.skipif(not stt.available(), reason="faster-whisper not installed (optional dep)")


def _silence_wav(seconds: float = 0.5, rate: int = 16_000) -> bytes:
    buf = io.BytesIO()
    with wave.open(buf, "wb") as w:
        w.setnchannels(1); w.setsampwidth(2); w.setframerate(rate)
        w.writeframes(b"\x00\x00" * int(rate * seconds))
    return buf.getvalue()


def _tone_wav(seconds: float = 0.5, rate: int = 16_000, hz: int = 440) -> bytes:
    frames = b"".join(struct.pack("<h", int(12000 * math.sin(2 * math.pi * hz * i / rate)))
                      for i in range(int(rate * seconds)))
    buf = io.BytesIO()
    with wave.open(buf, "wb") as w:
        w.setnchannels(1); w.setsampwidth(2); w.setframerate(rate); w.writeframes(frames)
    return buf.getvalue()


def test_empty_audio_is_empty_text_not_an_error():
    assert stt.transcribe(b"") == ""      # no model load, no exception


def test_oversized_recording_is_rejected_before_the_model_sees_it():
    with pytest.raises(stt.RecordingTooLarge, match="too large"):
        stt.transcribe(b"\x00" * (stt.MAX_BYTES + 1))


@needs_model
def test_undecodable_blob_is_NOT_reported_as_too_large():
    """PyAV raises ValueError subclasses on garbage. If the size guard also raised a plain
    ValueError, the endpoint answered 413 'recording too large' for a 9-byte junk upload - a wrong
    diagnosis that sends the user chasing a size limit. Keep the two failure modes distinguishable."""
    with pytest.raises(Exception) as ei:
        stt.transcribe(b"not audio")
    assert not isinstance(ei.value, stt.RecordingTooLarge)


def test_available_is_a_cheap_pure_gate():
    assert isinstance(stt.available(), bool)
    assert stt.available() is stt.available()    # lru_cached, no repeated import machinery


def test_warm_is_a_noop_when_the_dep_is_missing(monkeypatch):
    monkeypatch.setattr(stt, "available", lambda: False)
    stt.warm()      # must return immediately, never spawn a thread that explodes


@needs_model
def test_silence_transcribes_to_empty_string():
    assert stt.transcribe(_silence_wav()) == ""


@needs_model
def test_a_pure_tone_produces_no_hallucinated_words():
    # Whisper is notorious for inventing text on non-speech; vad_filter should suppress it.
    # If this ever fails loudly with a sentence, the VAD regressed - do not "fix" it by asserting the text.
    out = stt.transcribe(_tone_wav())
    assert len(out) < 30, f"hallucinated on a 440 Hz tone: {out!r}"


@needs_model
def test_browser_webm_opus_blob_decodes_in_process():
    """The real path: MediaRecorder posts webm/opus, PyAV decodes it, no ffmpeg binary involved."""
    av = pytest.importorskip("av")
    src = av.open(io.BytesIO(_tone_wav(seconds=0.4)), format="wav")
    buf = io.BytesIO()
    out = av.open(buf, mode="w", format="webm")
    ostream = out.add_stream("libopus", rate=48_000); ostream.layout = "mono"
    for frame in src.decode(audio=0):
        frame.pts = None
        for pkt in ostream.encode(frame):
            out.mux(pkt)
    for pkt in ostream.encode(None):
        out.mux(pkt)
    out.close(); src.close()
    blob = buf.getvalue()
    assert blob[:4] == b"\x1a\x45\xdf\xa3", "not an EBML/webm container"
    stt.transcribe(blob)      # must not raise: decoding is the assertion here
