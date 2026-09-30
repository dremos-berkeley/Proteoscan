#!/usr/bin/env bash
# Install ThermoRawFileParser 2.0.0-dev into the active conda env (or --prefix DIR).
#
#   conda activate peaksparity
#   bash scripts/install_trfp.sh              # self-contained .NET 8 build for this OS/CPU
#   bash scripts/install_trfp.sh --route dotnet   # framework-dependent build + conda dotnet-runtime 8
#
# Why not conda: bioconda's thermorawfileparser for macOS is 1.4.5 (noarch) and depends on
# mono, which conda-forge builds only for osx-64, so it cannot be installed natively on Apple
# Silicon. The upstream 2.0.0-dev release ships self-contained .NET 8 builds (no mono, no
# dotnet needed) for osx-arm64, osx (x86_64) and linux; on Linux bioconda already packages
# the same 2.0.0-dev version.
#
# macOS details handled here:
#   * the osx-arm64 apphost (ThermoRawFileParser) in the release zip is NOT code-signed, and
#     Apple Silicon kills unsigned arm64 executables ("Killed: 9"), so it is ad-hoc signed
#     with /usr/bin/codesign (part of macOS, no Xcode needed);
#   * the quarantine attribute is removed in case the zip was fetched by a browser.
set -euo pipefail

VERSION="2.0.0-dev"
BASE="https://github.com/CompOmics/ThermoRawFileParser/releases/download/v.${VERSION}"
# sha256 of the release assets as downloaded on 2026-09-30 (plain function, not an
# associative array: macOS ships bash 3.2)
sha_for() {
  case "$1" in
    osx-arm64) echo 6167aca5c4d44bbf35f986192bc249aba9760e138d13bba27791c2da591ac087 ;;
    osx)       echo 98bc806d5f45f32684d77d379c720bfd25aafe04cebcd4dc07f185fbc1942c59 ;;
    linux)     echo 19566762ce6759a93cee9aa4cef50de1e9ae2ad9078bd95826e8e733d4bb0d52 ;;
    net8)      echo c39221ef14a803cf59f45a46d6fe477ccf1fbf1007368364927dd56a84d01923 ;;
  esac
}

ROUTE="native"; PREFIX="${CONDA_PREFIX:-}"; VERIFY=1; FORCE=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --route) ROUTE="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --no-verify) VERIFY=0; shift ;;
    --force) FORCE=1; shift ;;
    -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown option $1" >&2; exit 2 ;;
  esac
done
[[ -n "$PREFIX" ]] || { echo "no conda env active; run 'conda activate peaksparity' or pass --prefix DIR" >&2; exit 1; }

OS="$(uname -s)"
if [[ "$ROUTE" == "dotnet" ]]; then
  ASSET="net8"
elif [[ "$OS" == "Darwin" ]]; then
  # choose by hardware, not by the (possibly Rosetta-translated) shell
  if [[ "$(sysctl -n hw.optional.arm64 2>/dev/null || echo 0)" == "1" ]]; then ASSET="osx-arm64"; else ASSET="osx"; fi
elif [[ "$OS" == "Linux" && "$(uname -m)" == "x86_64" ]]; then
  ASSET="linux"
else
  echo "no self-contained build for $OS $(uname -m); use --route dotnet" >&2; exit 1
fi

DEST="$PREFIX/share/thermorawfileparser-$VERSION"
WRAP="$PREFIX/bin/ThermoRawFileParser"
if [[ -e "$WRAP" && $FORCE -eq 0 ]]; then
  echo "$WRAP already exists ($("$WRAP" --version 2>/dev/null | tail -1 || echo 'not runnable')); use --force to replace" >&2
  exit 0
fi

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
ZIP="$TMP/trfp.zip"
echo "downloading ThermoRawFileParser-v.${VERSION}-${ASSET}.zip"
curl -fL --retry 3 -o "$ZIP" "$BASE/ThermoRawFileParser-v.${VERSION}-${ASSET}.zip"
if [[ $VERIFY -eq 1 ]]; then
  if command -v shasum >/dev/null; then GOT="$(shasum -a 256 "$ZIP" | cut -d' ' -f1)"; else GOT="$(sha256sum "$ZIP" | cut -d' ' -f1)"; fi
  if [[ "$GOT" != "$(sha_for "$ASSET")" ]]; then
    echo "sha256 mismatch for $ASSET: got $GOT, expected $(sha_for "$ASSET")." >&2
    echo "The release asset may have been re-uploaded; inspect it and rerun with --no-verify." >&2
    exit 1
  fi
fi
rm -rf "$DEST"; mkdir -p "$DEST" "$PREFIX/bin"
unzip -q "$ZIP" -d "$DEST"
# some zips have a top-level folder
if [[ ! -e "$DEST/ThermoRawFileParser.dll" ]]; then
  inner="$(dirname "$(find "$DEST" -name ThermoRawFileParser.dll | head -1)")"
  [[ -n "$inner" && "$inner" != "." ]] && { mv "$inner"/* "$DEST"/; }
fi

if [[ "$ASSET" == "net8" ]]; then
  command -v dotnet >/dev/null || { echo "dotnet not found: conda install -c conda-forge 'dotnet-runtime=8.0'" >&2; exit 1; }
  cat > "$WRAP" <<EOF
#!/usr/bin/env bash
exec dotnet "$DEST/ThermoRawFileParser.dll" "\$@"
EOF
else
  chmod +x "$DEST/ThermoRawFileParser"
  [[ -e "$DEST/createdump" ]] && chmod +x "$DEST/createdump"
  if [[ "$OS" == "Darwin" ]]; then
    xattr -dr com.apple.quarantine "$DEST" 2>/dev/null || true
    # ad-hoc sign every Mach-O that is not signed yet (the arm64 apphost ships unsigned)
    for f in "$DEST/ThermoRawFileParser" "$DEST/createdump" "$DEST"/*.dylib; do
      [[ -e "$f" ]] || continue
      codesign --verify "$f" 2>/dev/null || codesign --force --sign - "$f"
    done
  fi
  cat > "$WRAP" <<EOF
#!/usr/bin/env bash
exec "$DEST/ThermoRawFileParser" "\$@"
EOF
fi
chmod +x "$WRAP"
echo "installed: $WRAP -> $DEST"
"$WRAP" --version
