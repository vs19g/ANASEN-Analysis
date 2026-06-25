# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Git policy

- **Never push.** Only commit when explicitly asked. Remember only the last commit made.

## Build

ROOT macros (`.C`) compile via ACLiC inside ROOT — not linked into binaries. Standalone binaries (`EventBuilder`, `Mapper`) are built with `make` in `Armory/`.

```bash
cd Armory && make                              # build EventBuilder, Mapper
root -q -l -b -e '.L TrackRecon.C++O'         # pre-compile before parallel runs
```

Batch runs: `./run_17F.sh`, `./run_27Al.sh`, `./run_tr.sh`. No test suite — validation is visual via ROOT `TBrowser`.

## Runtime configuration

`TrackRecon.C::Begin()` reads all config from environment variables (`DATASET`, `reactiondata`, `source_vertex`, `Gain`, `CO2percent`, `A1C1_Z_SCALE`, `A1C1_Z_OFF_QQQ`, `A1C1_Z_OFF_SX3`, etc.). The batch scripts `export` these before launching ROOT. See `README.md` for the full pipeline and file reference.
