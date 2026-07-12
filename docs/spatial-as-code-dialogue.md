# Spatial as code - a speculative dialogue (English translation)

**Status:** vision / thinking-out-loud, NOT a grounded evidence base · **Scope:** where "biology as code" is most tractable (single-cell + spatial), what the language of "spatial as code" would be, what Heisenberg/Fisher bounds contribute, and what `rctd-py` could do that would be significant for the field.

> **Provenance.** This is a faithful English translation of a five-turn French dialogue (a back-and-forth exploring the framing). It is speculative by design - a source of framing and roadmap ideas, not verified findings. Nothing here is a claim to build on directly; treat the concrete `rctd-py` ideas at the end as hypotheses, not specs. For how much of it actually maps onto SpatialScribe, see the "Relevance to SpatialScribe" section at the bottom.

---

## Turn 1 - If we narrow to a sub-domain (single-cell RNA + spatial), could we do more?

Good instinct to narrow, because this is exactly a sub-domain where several of the locks on "biology as code" partly loosen. Let us think about what we could really do.

### Why this sub-domain is more favorable

Single-cell + spatial has three properties that bring it closer to an exploitable "as code" regime:

1. **It is massively quantified and already digital.** An scRNA-seq experiment is a cells x genes matrix. Spatial adds coordinates (x,y) - even (x,y,z) plus histology. We are in a dense, structured, high-dimensional data space: exactly the terrain of modern models. Unlike "design a cell," here we massively *measure* state rather than write a program to execute.

2. **The feedback loop is already shorter.** The major lock on "biology as code" was the design-build-test cycle at the scale of weeks. Here, much of the work is analytical: once the data are produced, iteration happens in silico, in minutes. The wet lab stays slow, but the "learn from the data" space is huge and fast.

3. **There is a readable notion of state.** A cell's gene expression is a reasonable proxy for its "program state" at time t. We do not see the execution, but we have snapshots - many snapshots.

### What we already know how to do (and which looks like "as code")

- **Cellular foundation models.** Geneformer, scGPT, scFoundation, Universal Cell Embedding: trained on tens of millions of cells, they learn a "grammar" of cell state. You can prompt them: mask genes and predict, or do zero-shot cell-type annotation. This is a learned abstraction layer over cell state - the analogue of the protein language model, but for the whole cell.
- **Perturbation prediction.** This is the most "as code" point: models (GEARS, scGen, CPA...) predict a cell's transcriptomic state after a never-tested knock-out or treatment. In other words: "if I edit this line (this gene), what does the behavior become?" - a compiler/simulator of the cellular response. Still imperfect, especially out of distribution, but exactly the missing brick.
- **Spatial reconstruction and communication inference.** Predict missing expression from histology, impute a full transcriptome from a limited panel (Tangram, gimVI), infer who talks to whom (ligand-receptor: CellChat, COMMOT). Spatial adds the dimension scRNA alone lacked: context, neighbors. Much of cell behavior is contextual - that is where spatial unlocks causality we did not have.

### Where "we could do more" - the real potential

Where the leverage is strongest in the near term:

1. **From static snapshots to dynamics.** Today we have photos. The frontier is reconstructing the trajectory - RNA velocity, spatial + time, lineage. The grail would be a model that, from a tissue, predicts its evolution: a tumor progressing, an embryo developing. We would move from a "debugger that inspects state" to an "emulator that runs the program." Playable because the data are numeric and spatio-temporal architectures exist.

2. **Tissue as a circuit, not the cell as an isolated unit.** Spatial lets you treat a tissue as a system: cells = nodes, proximity + signaling = edges. A graph neural network over a tissue models a "distributed program" running on a population. Conceptually much closer to "architecture as code" than "function as code" - you describe a spatial organization and its interaction rules.

3. **Design of spatially-targeted intervention.** Once you predict perturbations in context well, you can invert the problem: "what intervention brings this tissue from the diseased state to the healthy state?" - goal-conditioned generation. This is the true "intent -> compilation" of the domain.

### Where it still breaks (keep in mind)

- **Correlation != causation.** Snapshots give co-occurring states, not mechanisms. Perturbation helps, but the intervention space is gigantic and under-sampled.
- **Noise, dropouts, batch effects.** The signal is noisy, and current spatial trades off resolution (subcellular but few genes) vs coverage (whole transcriptome but multi-cell spot). We do not yet have "everything, everywhere, at cell resolution."
- **The arrow of time is hard to recover.** RNA velocity is an ingenious but fragile trick; inferring dynamics from statics stays ill-posed.
- **Out-of-distribution generalization.** A model trained on some tissues/species/conditions extrapolates poorly. The "compiler" is reliable only within its training domain.

### Summary

Single-cell + spatial is probably the sub-domain where "biology as code" is most advanced today, because it moves the problem from design (writing DNA, slow and costly) to analysis and prediction (learning from numeric matrices, fast and iterable). We already have the equivalent of learned abstraction layers (cellular foundation models) and partial simulators (perturbation prediction). The leap to "do much more" would come from three advances: dynamics (from snapshot to trajectory), tissue-as-system (the spatial GNN modeling a distributed program), and in-context intervention design.

Put differently: in single-cell alone, you read the state of a cell. Adding spatial, you read the distributed program of a tissue - the first biological object rich enough, numeric enough, and iterable enough for a real "compiler of the living" to start making sense at this scale.

*An angle to dig into: what would be the domain's own "fast verifier"? In software it is the compiler; here, is it a perturbation model reliable enough to filter in silico before the wet lab - and at what fidelity does that filter become economically transformative rather than just an approximate triage?*

---

## Turn 2 - Spatial as code: which Anthropic components would fit best?

A slightly counter-intuitive answer: the Anthropic components are not the "compiler of the living" - they are the orchestration and looping layer around that compiler. Claude does not replace scGPT, GEARS, squidpy, or an AlphaFold-like. Claude is the conductor that writes the pipeline, calls these specialized models as tools, reads their noisy outputs, and closes the generate -> verify -> correct loop. Keep this split in mind, or you pick the wrong brick.

Mapping, from most central to most peripheral.

### The core: the closed loop

- **Tool use (client tools) + Messages API.** The mechanism that is the "as code." Claude proposes a spatial spec or an intervention, calls the real verifier (a perturbation predictor, a spatial-neighborhood analysis, a statistical test), receives the `tool_result`, corrects. Claude responds with `stop_reason: "tool_use"`, your code executes the operation and returns a `tool_result`. Exactly the "design -> in-silico test -> learn" loop - except the verifier here is a learned bio model, not a formal compiler.
- **Code Execution (server Python/bash sandbox).** Probably the most underestimated component for this domain. Claude writes and directly executes scanpy / squidpy / anndata code in a sandbox: QC, clustering, spatial statistics, figures. The "code" of "spatial as code" is literally executed there.
- **Programmatic Tool Calling.** Decisive because of the data's nature. A single-cell + spatial matrix is heavy - you do not want it flowing through the context window. This mode lets Claude orchestrate tool calls inside the execution environment and manipulate large intermediates without loading them into context. Exactly what Claude for Excel does: read and modify sheets of thousands of rows without overloading the model's context. Swap "thousands of rows" for "hundreds of thousands of cells": same problem, same solution.

### The connective tissue: wiring the biological stack

- **MCP connector.** The answer to the "interoperability" argument, applied to bio. Rather than writing custom glue, wrap each specialized resource as an MCP server: a LIMS/Benchling, a spatial data warehouse, a compute cluster, and above all your internal perturbation model exposed as a tool. With MCP, the agent queries databases, integrates with APIs, and connects to services without writing bespoke tool implementations. The cellular foundation model then becomes "just" a tool Claude knows how to call.
- **Tool Search Tool.** A real bio stack is hundreds of methods (each imputation algo, each colocalization test, each model). You cannot preload them all. This tool lets Claude discover the right method on demand. Anthropic reports ~85% reduction in token use and a clear accuracy improvement on MCP evaluations when the tool library is large.

### Scale and autonomy

- **Batch API.** The design space is combinatorial: perturbations x cell types x spatial contexts. Batch serves mass in-silico screening - annotating millions of cells, filtering thousands of intervention candidates cheaply. You can even include MCP servers in Batch requests. This is what makes "cheap verifier at scale" economically credible.
- **Claude Agent SDK / Managed Agents.** To package the whole loop into an autonomous "design-build-test-learn" agent. The Client SDK makes you implement the tool loop yourself; the Agent SDK gives a Claude that runs the tools itself. Two key mechanisms: **hooks** (deterministic functions invoked at fixed points of the loop - your validation gates, the "CI/linter" equivalent of the pipeline: a candidate passes only if it satisfies a given toxicity/feasibility criterion), and **subagents** (a QC agent, a spatial-stats agent, an intervention-design agent - decomposition of the "distributed program").
- **Prompt caching + long context.** A spatial project's stable context is huge and reusable: protocols, panel definition, cell-type ontologies, markers, prior results, codebase. Cache all of it and vary only the query. This lowers the cost of each iteration of the slow loop.

### Model choice as "tiering"

An architect's reflex: Opus for hard reasoning (experimental design, hypothesis formulation, in-context intervention design); Sonnet for pipeline-generation throughput; Haiku for very-high-volume, low-cost triage (cell pre-annotation, coarse classification before fine filtering). Same principle as compute tiering in a real system.

### The honest recap

The winning architecture is not "Claude does the biology." It is Claude-as-conductor:
- it writes the spatial spec and the pipeline (Code Execution, Programmatic Tool Calling);
- it calls the real bio verifiers - your tissue GNNs, your perturbation predictors, your cellular foundation models - via tool use and MCP;
- it closes the loop and iterates (Agent SDK + hooks as validation gates);
- it scales the in-silico screen (Batch) before sending only the best candidates back to the wet lab.

Anthropic provides the orchestration, execution, connectivity, and looping layer. The "compiler/simulator of the living" layer - the model that actually predicts the transcriptome in spatial context - stays yours, and it is that layer which caps the fidelity of the whole system.

*The most "transformative" Anthropic component is no particular model - it is the tooled closed-loop mechanism (tool use + deterministic hooks). It only improves the system if the bio verifier you wire in is reliable enough to filter before the lab. The interesting thread: where to place the deterministic hooks - which gates must be hard rules (synthesis feasibility, toxicity) vs learned scores?*

---

## Turn 3 - What would the *language* of "spatial as code" be? Rust, TypeScript, a DSL?

Out of the plumbing. The real language-design question reveals that "spatial as code" cannot be Rust, TypeScript, or even a DSL in the classic sense. And the why is more interesting than the what.

### First: what is the program made of?

Rust and TypeScript share an invisible assumption: the machine is deterministic. Source -> binary -> reproducible execution. All their tooling (borrow checker, static types) rests on that. But the runtime of the living is stochastic, evolving, and rewrites your binary to save energy.

Radical consequence: a "spatial" program does not compile to a tissue, it compiles to a **distribution of possible tissues**. "Running" the program is sampling. The natural semantics is therefore not imperative - it is that of a probabilistic language (Stan, Pyro, Gen). You do not write "cell X does Y"; you write generative priors and constraints, and the observed tissue is a draw. This rules out the Rust/TS paradigm from the start: they describe instructions, not distributions.

### What we still steal from Rust and TypeScript

Each language has a central gift. The useful question is not "which to choose" but "which gift transposes."

**Rust's gift is static proof.** The borrow checker proves before execution that you do not violate memory. The transpose would be a "burden checker": linear types for the metabolic budget. Every construct costs ATP; the type system guarantees you do not exceed the cell's energy envelope - otherwise the runtime silences it. Ownership becomes metabolic ownership. Even better: annotate each construct with its evolutionary half-life - its resistance to being optimized away by selection. A type encoding "how many generations before the cell disables this gene because it is too costly." No language has thermodynamic or Darwinian types; here it is a necessity. But - the founding tension - you cannot prove biology statically. The checker is probabilistic: it returns confidence intervals, not theorems. Rust with a checker that doubts.

**TypeScript's gift is gradual typing over a chaotic substrate.** JS is a dynamic, untyped magma; TS adds contracts where it can guarantee them, and leaves `any` elsewhere. Probably the most honest analogy near-term. The cell is the JavaScript: living, disorderly, untyped reality. You annotate contracts only where a learned model can verify ("this region expresses this marker with this confidence"), and assume biological `any` everywhere else. Gradual typing over the living - progressively hardened, never total.

### The real paradigm is not linear, it is spatial

Rust and TS are line-by-line languages. A tissue has no linear execution order: everything happens in parallel, in space, as a field. The true ancestors are not general-purpose languages, they are:
- **shaders (GLSL):** you write one kernel per cell that reads the state of neighbors and morphogen fields. Space-parallel by nature.
- **HDLs (Verilog -> Cello):** you describe a desired spatio-temporal behavior, not steps; a synthesizer finds the implementation. Declarative, concurrent, spatial.
- **reaction-diffusion / cellular automata:** gradients, fields, local rules -> emergent global form.

So space and time are first-class primitives (coordinates, gradients, developmental trajectories), and the AST is not linear text - it is a **geometric graph**: cells = nodes, signaling = edges. The tissue-as-graph becomes the language's native data structure, not a library.

The contrast in two lines:

```
// Imperative (dead) - what we do NOT want
for cell in tissue: if neighbor.signal > k: cell.express(gene_A)

// Declarative field (living) - the intent + the constraint
field morphogen ~ diffuse(source=crypt_base, lambda=200um)
assert region(morphogen > theta): identity ~= Paneth +/- tol   @stability(>50 gen)
```

You do not code the loop. You declare a field and a target, and the "compiler" searches for the genetics that realizes it.

### The deepest break: the compiler is learned, so the language is differentiable

In software you debug by reading the error. Here the compiler is not formal, it is a learned world model. This suggests the language should be **differentiable by construction** (JAX-like semantics): you do not fix a bug, you descend a gradient from "desired tissue" toward "spec parameters." The compilation error becomes a loss. Debugging becomes optimization. "Type-checker" becomes "is the model confident this spec is realizable."

### The real twist: the source code is not text

Where the disruption is maximal. In any classic language, the language is primary and the compiler serves it. Here it is the **generative model that is primary**, and the language is only a human-readable interface to drive it.

The source of truth of "spatial as code" is probably not a text file. It is a coordinate in the latent space of a tissue world-model, plus a set of constraints. The "DSL" the human types is a lens - a query-and-constraint language projected onto that learned manifold, exactly as one queries embeddings. Text is the interface; the durable artifact is the latent point + the constraints that bound it.

This flips the whole stack:
- A cell's "type" is no longer a discrete class, it is a region in latent space with a tolerance. Subtyping becomes geometric proximity. Type-checking becomes "is this state inside the ball we asserted."
- The leaky modularity we deplored stops being a bug: on a continuous manifold there is no sharp boundary between "parts" - consistent with the living, where clean modules do not exist.
- Versioning becomes versioning of coordinates and constraints, not lines.

### The direct answer

Not Rust, not TypeScript, not a classic textual DSL. The language of "spatial as code" is a **constraint-and-intent lens laid over a differentiable world-model** - probabilistic in semantics (you compile to distributions), gradual in typing (TS-like, over a living substrate), spatial and field-native in structure (shader/HDL-like), endowed with resource types for metabolic cost and evolutionary stability (Rust's gift, but with a checker that doubts), and whose real "source code" is a latent point, not text.

The pretty paradox: the better the underlying model, the less expressive the language needs to be - it slides toward mere specification of objectives and guardrails, because the model carries all the semantics. At the limit, "spatial as code" is no longer programming, it is constrained steering of a world-model, and the DSL is just the surface where you write the invariants you refuse to see violated.

And this reframes where a system like Claude stands: it is neither the language nor the world-model - it is the bidirectional translator between human intent and latent steering. It takes "I want an intestinal crypt that resists inflammation," writes it into the constraint lens, reads what the model answers, and negotiates the trade-offs (metabolic budget vs stability vs function). The compiler carries the biology; the language carries the invariants; the agent carries the translation and the loop.

*Next thread: in a language whose verifier doubts, which control primitive do we refuse to make probabilistic? You necessarily need a handful of hard, non-negotiable invariants (biocontainment, kill-switch, non-transmissibility) that stay proofs and not scores - otherwise you steer a world you no longer bound. Where that boundary between learned-declarative and hard-proved falls is perhaps the design decision of the whole domain.*

---

## Turn 4 - Can Heisenberg's ideas be useful here?

Yes - but be surgical, because Heisenberg is the physicist most often cited badly. The temptation is to slap "everything is uncertain in biology, therefore uncertainty principle." That is pop-physics. The Heisenberg ideas that actually bite here are more precise. Three are structuring, one is a trap, and one design consequence is genuinely new.

### The real conjugate pair: state vs trajectory

Heisenberg does not say "everything is blurry." He says something stronger: certain pairs of quantities are non-co-measurable - knowing one precisely forbids knowing the other. Position and momentum. So the useful question is not "is there noise?" but "what is the conjugate pair of single-cell?"

And it exists, physically, not by analogy: the **state and the trajectory** of the same cell. scRNA-seq destroys the cell to read its transcriptome. You can know a cell's full state (at the price of its death) or watch it become what it becomes (at the price of not lysing it) - never both on the same individual. Exact position, or dynamics: not both.

This reframes something touched earlier. RNA velocity was called "an ingenious but fragile trick." Heisenberg explains why it is *structurally* fragile, not just technically: it tries to infer the conjugate variable (momentum, dynamics) from a single position snapshot. Exactly the operation the principle forbids doing cleanly. Velocity is not a bad algorithm - it is an attempt to bypass a non-commutation. So the "spatial as code" language must treat "current state" and "dynamics" as a conjugate pair sampled separately over populations, never jointly on an individual. A semantic constraint, not an implementation one.

### The measurement effect as a tax on the loop

The thread's through-line is the generate -> verify -> correct loop and the dream of a fast verifier. Heisenberg taxes it: to measure is to perturb. The wet lab that anchors the in-silico model is destructive and non-repeatable - you never measure the same cell twice. The verifier does not return "the truth about this individual," it returns "the truth about a sacrificed individual, in a supposedly equivalent population."

Concrete architectural consequence: the loop can never be longitudinal at the unit scale. It is statistical by construction - you infer a population's dynamics from cohorts destroyed at different times. This changes what "correct" means: you do not correct an observed trajectory, you update a distribution from destructive slices. The deterministic hook we wanted in the Agent SDK must validate on populations, never on tracked cells.

### Matrix mechanics: the license for anti-realism

The deepest point, and it retroactively legitimizes the previous turn. Heisenberg's real move in 1925 was not uncertainty - it was matrix mechanics: stop describing the orbits of unobservable electrons, describe only the relations between observables (spectral lines, transition amplitudes). An epistemology: speak only of what you measure, stay silent on the hidden essence.

That is exactly the status of the learned world-model. A cellular foundation model does not model "the true cell"; it models relations between observed states. "The type of a cell is no longer a discrete class but a region in latent space with a tolerance" was Heisenberg without knowing it: you define an entity by its observable relations, not by a hidden mechanistic identity. The leaky modularity we deplored is no longer a defect - it is the honest admission that there are only observables and their transitions.

The corollary, his most famous line: "what we observe is not nature itself, but nature exposed to our method of questioning." Translation for "biology as code": your latent space does not have the shape of nature, it has the shape of the assay. Batch effects, dropout, platform artifacts - the "compiler" baked them in. Not noise to clean, but constitutive. So the language needs a primitive no classic language has: **the provenance of the questioning method as part of the type**. A constraint is not "this region is Paneth," it is "this region is Paneth as exposed by MERFISH at this depth." The `@stability(>50 gen)` should be flanked by an `@observed_via(...)`.

### The trap: there is no biological h-bar

Where you must be honest, or become a charlatan. Heisenberg's uncertainty is ontological - fundamental, not an instrument defect. No better camera beats dx*dp >= h/2.

Most biological "uncertainties" are epistemic/instrumental: the resolution-vs-coverage trade-off (subcellular or whole-transcriptome) could be beaten by better tech. Treating it as a law of nature would be the classic error. The state/trajectory pair is closer to fundamental - but even it yields if non-destructive live imaging advances. So Heisenberg is an intuition pump here, not a theorem to import. There is no Planck constant of the living.

What matters is that the rigorous theorem exists anyway - it is just not Heisenberg's. The honest mathematical form of "biological uncertainty" is information-theoretic: the **Cramer-Rao bounds** (the Fisher information a finite number of cells can carry about a parameter), and **rate-distortion theory** (how many genes x how much spatial resolution a channel can transmit). That is the real conjugate bound, and it is fundamental once the protocol is fixed. Heisenberg gives the intuition "there are pairs you cannot tighten together"; Fisher/Shannon give the number. The good design reflex: for each conjugate pair of the domain, seek its Cramer-Rao bound rather than invoking an imaginary h-bar.

### The design consequence: a conjugate-aware type system

All of this converges on a new language primitive. In Rust, a value has a type and a value. Here, every measured quantity carries a value and a resolution, and the type-checker must forbid asserting two conjugate quantities at a joint precision exceeding what is co-measurable.

```
// forbidden by the checker: joint over-precision of a conjugate pair
assert cell.state == PanethMature @res(single_cell)
   and cell.fate  == Absorptive    @res(single_cell)   // X non-co-measurable

// allowed: the conjugate is only asserted at population level
assert cell.state == PanethMature @res(single_cell)
   and population.fate_distribution ~ ... @res(cohort)  // OK
```

A verifier that catches neither memory errors (Rust) nor type errors (TS), but **over-claiming errors**: assertions that pretend to jointly know things the measurement channel cannot deliver together. The metabolic "burden checker" of the previous turn gains a twin: an epistemic conjugate checker. It does not say "this program will crash," it says "this program claims to know more than is measurable." Probably Heisenberg's most disruptive contribution to this language: a type system whose errors are epistemological.

### Back to the hard-invariant boundary

The previous turn ended on: where does the line between learned-declarative and hard-proved (biocontainment, kill-switch) fall? Heisenberg settles part of the debate: you can only erect as a hard invariant what you can verify without destroying the system. A kill-switch whose only verification is destructive is not a guarantee - it is an a-posteriori belief about already-dead individuals. Non-co-measurability bounds the set of invariants that can be hard. A safety invariant is truly hard-provable only if it falls in the subspace of quantities observable continuously, in vivo, without perturbing. So the declarative/proved boundary is not only a design choice - it is partly dictated by the conjugate structure of the observables.

*Next thread: this suggests a credible safety requirement in "biology as code" forces you to design the organism to be non-destructively observable - to make continuous auditability an upstream design constraint, not a downstream check. Safety is not inspected, it is compiled into the living. A fairly radical inversion of how we think about software verification - and there, curiously, it is Bohr and complementarity rather than Heisenberg who take over.*

---

## Turn 5 - With all this, and taking the `rctd-py` codebase, what significant thing could we do for the spatial field?

An excellent anchor - because RCTD is not a random repo: it is precisely the tool that lives on the fault line the previous turns spent describing. And there is a lever in this repo that, I believe, no one has pulled yet.

### What this repo really is (and the hidden lever)

RCTD deconvolves: a spatial spot contains a mix of several cells, and the algorithm recovers the cell-type proportions by leaning on an scRNA reference. In other words - it is the **un-mixing operation of the resolution/coverage conjugate** named in the Heisenberg turn. The spot gives you whole-transcriptome coverage at the price of single-cell resolution; RCTD tries to recover the lost composition. It is the computational form of the domain's fundamental trade-off. And its platform-effect correction step (the sigma, the `fit_bulk`) is literally Heisenberg-matrix-mechanics: the model does not pretend to see "the true expression," it models the relation between two questioning methods (the reference and the spatial assay) that distort genes differently. The README says it itself: the 0.27% deviation from R comes entirely from the platform-effect estimation, not the solver. The assay-shaped latent, in one benchmark sentence.

But here is the lever. The original R version is not differentiable. By porting it to PyTorch, the author - perhaps without aiming for it - made the whole deconvolution **differentiable, GPU-native, and validated bit-identical at whole-transcriptome scale**. Not a speedup. A change of nature. Recall the thesis two turns earlier: "the language of spatial as code is differentiable by construction." Here is a real, published, exact object that just crossed into that category. Everything below exploits that crossing.

### Flagship idea - a deconvolution that refuses to over-claim

RCTD's deep problem, the one no speed fixes: deconvolution is an under-determined inverse problem. A pixel's counts are compatible with a *distribution* of compositions, not a single one. Yet RCTD crushes that distribution into a point (NxK weights) plus a hard singlet/doublet/reject classification. The kicker: the README *detects* the ambiguity - pixels disagreeing with R have a margin < 0.05 between the two best types - and it throws it away as a binary reject. The most precious information of the computation is produced then erased.

The fix is almost free once you are in PyTorch. The per-pixel solver is a constrained IRWLS - a maximum likelihood that already forms the weighted curvature. The Hessian of the negative log-likelihood at the MLE is the observed Fisher information; its inverse is the asymptotic covariance of the composition estimate. So with autograd you get, per pixel, not "type A at 0.6, type B at 0.4" but a **calibrated posterior over composition**, and a **Cramer-Rao bound** on the precision any method could reach given that pixel's UMI depth and the separability of the types in that reference.

Exactly the payoff of the Heisenberg turn - not an imaginary biological h-bar, but the real information bound (Fisher/Cramer-Rao) that is the honest form of "the uncertainty of the living." And it is the conjugate checker of the type system, instantiated in a running tool: RCTD would emit, per pixel, "these two types are information-theoretically non-separable here - do not believe a confident call." The semantics no classic language has (epistemological errors, not memory errors) becomes a column in the output h5ad.

The honest subtlety, and where the real science is: the weights live on a simplex (positive, sum to 1). The Gaussian Laplace approximation breaks at the boundary (when a weight tends to 0) - and that boundary is precisely the singlet/doublet frontier. So you need a posterior adapted to the constraint: logistic-normal or Dirichlet on the simplex, or a projected covariance. Not an annoying detail but the heart of the matter - the "simplex interior = confident mixture" vs "boundary = ambiguous identity" distinction is geometrically the distinction RCTD already hard-encodes, made continuous and calibrated.

Why it matters for the field, not just elegant: it changes the *type* of RCTD's output, from "an answer" to "a calibrated belief with a resolution limit." All downstream analysis - differential composition between conditions, niche detection, enrichment tests - currently treats the weights as truth. Giving them propagatable error bars is of systemic importance, and almost nobody delivers it at whole-transcriptome scale on GPU. `rctd-py` is the only substrate where it is within a Hessian's reach.

### Idea 2 - joint deconvolution: tissue-as-graph in the loss

RCTD deconvolves each pixel independently. It ignores that neighboring pixels are spatially correlated - a hepatocyte is probably surrounded by hepatocytes. A gaping hole given everything said about "tissue as distributed program." And crucially: the neighborhood is the missing information that makes the inverse under-determined. A mixed spot alone is ambiguous; its neighbors provide the constraint that lifts the ambiguity.

Because it is differentiable, you simply add a spatial-coupling term to the objective - a graph Laplacian (smoothing neighbor compositions) or, more ambitious, a learned prior via message-passing on the tissue graph - and backpropagate jointly. Per-pixel deconvolution becomes joint tissue deconvolution. The "GNN on the tissue graph" primitive made concrete and useful: the neighborhood prior refines exactly the pixels Idea 1 flags as ambiguous. The two ideas compose - the uncertainty tells you where the neighborhood prior should weigh most, and the prior resolves what the isolated pixel could not.

### Idea 3 - designing the next measurement: information-driven panel optimization

The destructive-measurement thread that closed the Heisenberg turn. Targeted platforms (Xenium, MERFISH) impose a gene budget - typically 300 to 500. Which to choose? Today, by hand or by marker heuristics.

With the per-pixel Fisher information of Idea 1, you compute, for a candidate gene set, the expected reduction in composition uncertainty for the cell types you care about. Optimal experimental design - maximize Fisher information, minimize Cramer-Rao - posed over gene panels. Rate-distortion applied to panel choice. Concretely: `rctd-py` + Fisher becomes a panel recommender that answers "to resolve TH1 vs TH2 in this tissue, these 40 genes are information-optimal."

The idea posed at the very end - "safety/observability compiles upstream" - becomes a tool: you no longer just analyze the assay you ran, you design the assay so the types you target become resolvable. This closes the loop from analysis back to measurement design - what the destructive-measurement constraint forces, since you will never remeasure the same cell. Very high practical value, and real money behind panel choice.

### The synthesis: from tool to recognition network

Stack the three and `rctd-py` stops being "RCTD but faster" and becomes a calibrated, spatially-coupled, measurement-design-aware inference engine. But the long horizon is prettier, and it loops back to the differentiable world-model of two turns earlier.

A probabilistic RCTD is the **recognition network** - the encoder - of a generative tissue model. The natural next step: make the reference itself learnable. Since everything is differentiable, you can backpropagate spatial fit quality into the reference profiles, under a prior - and thus *learn* the alignment between the two assay latent spaces instead of correcting it with a mere multiplicative factor. At that point `rctd-py` becomes the decoder-half of the "differentiable tissue world-model," its first concrete brick, built on a bit-identical-validated foundation rather than sand. And the substrate is already ready for the agentic loop: JSON output "for AI agents," CLAUDE.md, structured CLI. The "Claude-as-conductor" of the early turns already has its handles, nothing to wire.

### The honest limit - and the deepest trap

Two real fragilities. The Fisher/Laplace posteriors are asymptotic and fragile at the simplex boundary; the neighborhood prior can over-smooth and erase the true sharp frontiers - that is, the very tissue architecture you care about. So the coupling must be edge-aware, which loops back to needing the uncertainty to decide where to trust the smoothing. Everything holds together, but everything depends on one link.

And the trap, the deepest of all: **a well-calibrated posterior under a misspecified reference is worse than an honest point estimate - because it looks trustworthy.** Calibration is not correction. Their own benchmark contains the perfect example: the cervix deconvolved against a normal atlas for lack of a cancer atlas. A wrong reference produces confident, wrong posteriors.

Which gives the next thread, and it closes the arc of the whole conversation. The real frontier is not "deconvolve better." It is: **can the spatial data tell you the reference is misspecified** - a calibrated "the questioning method does not fit this tissue" signal - instead of silently producing confident nonsense? That is the honest end of Heisenberg's reasoning: knowing when your questioning method has distorted things beyond what the model captures. A misspecified generative model leaves a fingerprint in the residuals; a differentiable, probabilistic deconvolution could surface it. A deconvolution that can say "my reference is wrong for this tissue" - that would be the truly significant contribution. And it is precisely because the repo became differentiable that it is, for the first time, within reach.

---

## Relevance to SpatialScribe (the AI-annotation decision platform)

Honest triage of the above against what SpatialScribe actually is - a self-serve annotation copilot whose differentiator is a per-cell **confidence / abstention layer** (`Unassigned` / `Ambiguous` / `Unresolvable` / `Uncertain` / `Novel`), a **panel-adequacy check** (H3), and the thesis **reference choice > method choice**.

### Directly relevant - the one load-bearing nugget

- **rctd-py + Fisher information -> calibrated per-cell composition uncertainty (Idea 1) is the theoretical grounding for SpatialScribe's abstention layer.** The five abstention classes are today driven by heuristics/thresholds ([`annotation_qc_thresholds.yaml`](research/annotation_qc_thresholds.yaml)). Idea 1 gives them a principled backbone: a Cramer-Rao "these two types are information-non-separable here" signal maps almost one-to-one onto `Ambiguous` / `Unresolvable` / `Uncertain`. This is the most transferable idea in the whole dialogue, and it lands exactly on the platform's differentiator. SpatialScribe already routes annotation through RCTD (`rctd-annotation` / `rctd-py` skills, `annotate.py` spec) - so it connects to code you already run.
- **Reference-misspecification detection (the Turn-5 "deepest trap") is a concrete, high-value feature idea.** It is the exact restatement of SpatialScribe's own "reference choice > method choice" thesis and its two baked-in gotchas: the 5K-panel melanocyte gap and the "never let Claude invent numbers / gene presence != detectability" rule. A calibrated "my reference is wrong for this tissue" signal is precisely the abstention the platform wants to surface to a wet-lab user who supplies a mismatched atlas. Worth an entry in [`EXTENSIONS.md`](EXTENSIONS.md).
- **Panel-adequacy as *design*, not just diagnosis (Idea 3).** SpatialScribe's H3 currently *diagnoses* an inadequate panel (the melanocyte gap). Idea 3 inverts it into a Fisher-information panel *recommender*. That is a natural stretch/roadmap extension of an existing feature, not a new pillar.

### Validates the current architecture (no new action)

- **Turn 2's Anthropic-components mapping is essentially SpatialScribe's architecture already:** Agent SDK tool schema wrapping `analysis/` functions (whitelisted, not arbitrary exec), the action log as the deterministic record, hooks as validation gates, model tiering. The dialogue confirms the design choices rather than adding to them. Reassuring, not actionable.

### Inspiring framing, but out of scope for a hackathon annotation tool

- The grand "biology as code" arc - compiling DNA, intervention design, kill-switch/biocontainment, the differentiable tissue world-model, the "spatial as code" DSL, latent-point-as-source - is a research moonshot, not a self-serve copilot feature. Keep it as narrative/vision (a good "why this matters" for a pitch or README intro), but do not let it pull MVP scope. SpatialScribe reads a tissue's state and flags when it cannot trust its own labels; it is not steering a world-model.

### Bottom line

~80% of the dialogue is stimulating vision that should stay as framing. The remaining ~20% is unusually well-aimed: the rctd-py uncertainty-quantification and reference-misspecification ideas are the theoretical spine of the exact thing that makes SpatialScribe more than a wrapper - its confidence/abstention layer. If any single idea here becomes a roadmap item, it is **"calibrated per-cell composition uncertainty + a reference-fit signal, surfaced as abstention."**
