# prideQC OpenMS develop container fix v2

This supersedes the previous source-build container overlay.

Fixes:
- OpenMS source defaults to `https://github.com/OpenMS/OpenMS.git` / `develop`.
- `build_docker.sh` resolves the requested branch/tag to a full 40-character SHA before Docker starts.
- The Docker `pyopenms-builder` now fetches that immutable SHA (`OPENMS_GIT_REF`), not the moving requested branch (`OPENMS_GIT_REQUESTED_REF`). This removes the race where `develop` advances between `git ls-remote` and the Docker fetch.
- `nanobind-backend==1.0.0` is reinstalled together with the local pyOpenMS wheel after each `uv sync` that can prune it.
- build scripts remain executable.

The requested ref is still recorded in `/opt/prideqc/openms-build-info.txt` for provenance; the resolved SHA is the source actually built.
