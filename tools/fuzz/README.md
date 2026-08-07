# Differential fuzzing

picard-rs fuzzes tools through their arguments, because it ships tools. This repository ships a
library, so what there is to fuzz is the parsers a hostile file reaches first: `ITF8`, `LTF8` and
the CRAI index line.

Each is a **pure function of its bytes**. No state, no environment, nothing to seed. That is what
makes a divergence between the two sides a bug rather than a difference of setup, and it is why
this fuzzer needs neither a warm JVM nor coverage instrumentation to be worth running.

## Run

```sh
python3 tools/fuzz/run_fuzz.py --parser itf8 --iterations 300 --seed 1
```

The two sides are `FuzzDriver.java`, run in the pinned container, and
`crates/htsjdk-cram/examples/differential_fuzz.rs`. Both print one line per input, so the
comparison is a diff rather than a report.

## What a finding is

`findings/` holds JSON, written only when the two sides disagree: the shortest input that still
shows the divergence, what each side did with it, and the input it was minimised from. Bytes are
dropped from the end for as long as the divergence survives.

Findings are **evidence, not goldens**. One produced anywhere other than real x86-64 CI cannot be
quoted (decision 0008), which is why the job publishes them as an artefact rather than committing
them.

## The seed is an argument

Not the clock. A finding is reproduced by rerunning the command that found it, and a run that finds
nothing says so about a corpus anyone can regenerate.
