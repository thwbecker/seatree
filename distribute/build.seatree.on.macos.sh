#!/usr/bin/env bash

# SEATREE macOS Builder/Packager
# - Stages a portable folder and creates a DMG for distribution
# - Does NOT modify source code; assumes rpaths are already correct
#
# Requirements (preinstalled by user):
#   - Xcode Command Line Tools
#   - Homebrew packages: gtk4, gawk, ghostscript
#   - Python 3 with pygobject and matplotlib
#
# Usage:
#   ./build.seatree.on.macos.sh               # builds DMG to dist/Seatree-<date>.dmg
#   ./build.seatree.on.macos.sh --name Seatree-1.0 --out dist
#   ./build.seatree.on.macos.sh --no-dmg       # stage folder only (no DMG)

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
OUT_DIR="${SCRIPT_DIR}/dist"
VOL_NAME="SeatreeApp"
APP_NAME="Seatree"
DMG_NAME="Seatree-$(date +%Y%m%d)"
MAKE_JOBS="$(sysctl -n hw.ncpu 2>/dev/null || echo 4)"
CREATE_DMG=1
STAGE_DIR="${OUT_DIR}/${APP_NAME}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --name)
      shift; DMG_NAME="${1:-${DMG_NAME}}";;
    --out)
      shift; OUT_DIR="${1:-${OUT_DIR}}"; STAGE_DIR="${OUT_DIR}/${APP_NAME}";;
    --no-dmg)
      CREATE_DMG=0;;
    -h|--help)
      echo "Usage: $0 [--name <dmg_basename>] [--out <output_dir>] [--no-dmg]"; exit 0;;
    *)
      echo "Unknown option: $1"; exit 1;;
  esac
  shift || true
done

echo "==> SEATREE macOS build"
echo "ROOT_DIR: $ROOT_DIR"
echo "OUT_DIR:  $OUT_DIR"
echo "STAGE:     $STAGE_DIR"
echo "DMG base:  $DMG_NAME"

# 1) Sanity checks (no changes are made)
echo "==> Checking prerequisites"
missing=()
command -v hdiutil >/dev/null || missing+=("hdiutil")
command -v python3  >/dev/null || missing+=("python3")
command -v gawk     >/dev/null || missing+=("gawk (brew install gawk)")
command -v gs       >/dev/null || missing+=("ghostscript 'gs' (brew install ghostscript)")
pkg-config gtk4 --exists 2>/dev/null || missing+=("gtk4 runtime (brew install gtk4)")

if [[ ${#missing[@]} -gt 0 ]]; then
  echo "ERROR: missing requirements:" >&2
  for m in "${missing[@]}"; do echo "  - $m" >&2; done
  echo "Install missing tools, then re-run." >&2
  exit 1
fi

# Optional Python deps check (non-fatal):
python3 - <<'PY' || true
import sys
missing = []
for mod in ("gi", "matplotlib"):
    try:
        __import__(mod)
    except Exception:
        missing.append(mod)
if missing:
    print("WARNING: Missing Python packages:", ", ".join(missing))
    print("Install via: pip3 install matplotlib pygobject (or brew + pip setup)")
PY

# 2) Rebuild HC/GMT consumers if desired (uses existing Makefiles/rpaths)
echo "==> Building HC (skipped - use pre-built binaries)"
# Disabled due to makefile path issues - assumes HC binaries are already built
# if [[ -f "${ROOT_DIR}/modules/mc/hc/Makefile" ]]; then
#   (cd "${ROOT_DIR}/modules/mc/hc" && make -j"${MAKE_JOBS}" clean all)
# else
#   echo "Skip: modules/mc/hc/Makefile not found"
# fi

# 3) Stage payload
echo "==> Staging payload"
rm -rf "${STAGE_DIR}"
mkdir -p "${STAGE_DIR}"

# Copy project files needed to run locally from any mount point
RSYNC_EXCLUDES=(
  --exclude ".git" --exclude ".DS_Store" --exclude "*.o" --exclude "*.pyc" --exclude "__pycache__"
  --exclude "objects" --exclude "*.AppDir" --exclude "distribute" --exclude "*.dSYM" --exclude "dist"
  --exclude "tmp"
)
rsync -a --delete "${RSYNC_EXCLUDES[@]}" \
  "${ROOT_DIR}/python3" \
  "${ROOT_DIR}/modules" \
  "${ROOT_DIR}/seatree" \
  "${ROOT_DIR}/"*.md \
  "${STAGE_DIR}/"

# Copy optional files if they exist
for f in gmt-4.5.18 netcdf-c-4.9.3-rc1 COPYING COPYRIGHT README; do
  [[ -e "${ROOT_DIR}/$f" ]] && rsync -a "${ROOT_DIR}/$f" "${STAGE_DIR}/"
done

# Add a double-clickable launcher for Finder users
cat > "${STAGE_DIR}/Launch Seatree.command" <<'SH'
#!/usr/bin/env bash
DIR="$(cd "$(dirname "$0")" && pwd)"
# Provide a stable root for portable configs
export APPIMAGE_ROOT="$DIR"
# Ensure Homebrew user tools are available if needed
export PATH="/opt/homebrew/bin:/usr/local/bin:$PATH"
cd "$DIR"
exec ./seatree
SH
chmod +x "${STAGE_DIR}/Launch Seatree.command"

# 3b) Rewrite configuration files to be portable (relative paths)
echo "==> Rewriting configuration files for portability"

# Base paths inside staged folder
PY_BASE="python3"    # SEATREE base (where conf/ resides)
CONF_DIR="${STAGE_DIR}/${PY_BASE}/conf"

# Ensure conf subdirs exist
mkdir -p "${CONF_DIR}/hc" "${CONF_DIR}/conman" "${CONF_DIR}/syn2d"

# hcConf.xml: use APPIMAGE_ROOT for portable absolute resolution
cat > "${CONF_DIR}/hc/hcConf.xml" << 'XML'
<?xml version="1.0" encoding="UTF-8"?>
    <FlowConfiguration>
     <hcPath>APPIMAGE_ROOT/modules/mc/hc/bin</hcPath>
    </FlowConfiguration>
XML

# conmanConf.xml: use APPIMAGE_ROOT for portable absolute resolution
cat > "${CONF_DIR}/conman/conmanConf.xml" << 'XML'
<?xml version="1.0" encoding="UTF-8"?>
    <ConManConfiguration>
     <conmanPath>APPIMAGE_ROOT/modules/mc/ConMan/src</conmanPath>
    </ConManConfiguration>
XML

# syn2dConf.xml: use APPIMAGE_ROOT for portable absolute resolution
cat > "${CONF_DIR}/syn2d/syn2dConf.xml" << 'XML'
<?xml version="1.0" encoding="UTF-8"?>
    <Syn2DConfiguration>
     <chkbdPath>APPIMAGE_ROOT/modules/seismo/syn2d/makemodel/bin</chkbdPath>
     <makedataBinPath>APPIMAGE_ROOT/modules/seismo/syn2d/makedata/bin</makedataBinPath>
     <invertBinPath>APPIMAGE_ROOT/modules/seismo/syn2d/inversion/bin</invertBinPath>
    </Syn2DConfiguration>
XML

# conf.xml: set GMT path relative and keep module imports; directories relative to staged root
cat > "${CONF_DIR}/conf.xml" << 'XML'
<?xml version="1.0" encoding="UTF-8"?>
    <SEATREEConfiguration>
     <modules>
      <module>
       <importName>seatree.modules.hc.flowCalc</importName>
       <className>FlowCalc</className>
       <directory>python3/seatree/modules/hc</directory>
      </module>
      <module>
       <importName>seatree.modules.syn2d.syn2d</importName>
       <className>Syn2D</className>
       <directory>python3/seatree/modules/syn2d</directory>
      </module>
      <module>
       <importName>seatree.modules.conman.conman</importName>
       <className>ConMan</className>
       <directory>python3/seatree/modules/conman</directory>
      </module>
      <module>
       <importName>seatree.modules.larry.invert</importName>
       <className>Invert</className>
       <directory>python3/seatree/modules/larry</directory>
      </module>
      <module>
       <importName>seatree.modules.nonLinLoc.NonLinLoc</importName>
       <className>NonLinLoc</className>
       <directory>python3/seatree/modules/nonLinLoc</directory>
      </module>
      <module>
       <importName>seatree.modules.larry3d.larry3d</importName>
       <className>larry3d</className>
       <directory>python3/seatree/modules/larry3d</directory>
      </module>
      <module>
       <importName>seatree.modules.eqdyna.eqdyna</importName>
       <className>EQDYNA</className>
       <directory>python3/seatree/modules/eqdyna</directory>
      </module>
     </modules>
     <gmtPath></gmtPath>
     <convertPath>/usr/bin</convertPath>
    </SEATREEConfiguration>
XML

echo "==> Staged at: ${STAGE_DIR}"

# 4) Build a self-contained Seatree.app bundle inside dist (no source changes)
APP_DIR="${OUT_DIR}/Seatree.app"
APP_CNT="${APP_DIR}/Contents"
APP_MAC="${APP_CNT}/MacOS"
APP_RES="${APP_CNT}/Resources"
APP_FWK="${APP_CNT}/Frameworks"

echo "==> Creating app bundle at: ${APP_DIR}"
rm -rf "${APP_DIR}"
mkdir -p "${APP_MAC}" "${APP_RES}" "${APP_FWK}"

# Copy staged Seatree tree into app Resources
rsync -a --delete "${STAGE_DIR}/" "${APP_RES}/Seatree/"

# Copy icon if it exists
if [[ -f "${ROOT_DIR}/Seatree.icns" ]]; then
  cp "${ROOT_DIR}/Seatree.icns" "${APP_RES}/"
  echo "  - copied app icon"
fi

# Vendor Homebrew runtimes inside app (GTK, GI, Python, gawk, ghostscript)
echo "==> Vendoring Homebrew runtimes into app"
BREW_PREFIX="$(brew --prefix)"
if [[ -z "${BREW_PREFIX}" ]]; then
  echo "ERROR: Homebrew not found. Cannot vendor runtimes." >&2
  exit 1
fi

vendor_dir="${APP_RES}/vendor/opt"
mkdir -p "${vendor_dir}"

FORMULAE=(gtk4 glib pango cairo harfbuzz fribidi gdk-pixbuf gobject-introspection gawk ghostscript python@3.13 netcdf netcdf-fortran hdf5 jpeg libpng libtiff zlib gettext gmt)
for f in "${FORMULAE[@]}"; do
  src="${BREW_PREFIX}/opt/${f}"
  if [[ -e "$src" ]]; then
    echo "  - copying ${f}"
    rsync -aL "$src/" "${vendor_dir}/${f}/" || echo "WARNING: Failed to copy ${f}, continuing..."
  else
    echo "WARNING: ${f} not found at ${src}; skipping"
  fi
done

# Python packages: copy key site-packages from builder's Python AND system gi/cairo
SITE_OUT="${APP_RES}/python-site"
export SITE_OUT
mkdir -p "${SITE_OUT}"

# Copy gi and cairo from Homebrew's system Python 3.13 site-packages
# Use -L to follow symlinks and copy actual files
BREW_PY_SITE="${BREW_PREFIX}/lib/python3.13/site-packages"
if [[ -d "${BREW_PY_SITE}/gi" ]]; then
  echo "  - copying gi from ${BREW_PY_SITE}/gi (following symlinks)"
  rsync -aL "${BREW_PY_SITE}/gi" "${SITE_OUT}/"
fi
if [[ -d "${BREW_PY_SITE}/cairo" ]]; then
  echo "  - copying cairo from ${BREW_PY_SITE}/cairo (following symlinks)"
  rsync -aL "${BREW_PY_SITE}/cairo" "${SITE_OUT}/"
fi

# Copy other Python packages needed
python3 - <<'PY' || true
import sys, site, shutil, os
print('Python:', sys.version)
paths = []
try:
    paths = site.getsitepackages() + [site.getusersitepackages()]
except Exception:
    paths = [p for p in sys.path if p.endswith('site-packages')]
seen = set()
targets = ['PIL','Pillow','matplotlib','numpy']
out = os.environ['SITE_OUT']
for base in paths:
    if not os.path.isdir(base):
        continue
    for t in targets:
        src = os.path.join(base, t)
        if os.path.isdir(src) and t not in seen:
            dst = os.path.join(out, t)
            print('  - copy', src, '->', dst)
            shutil.copytree(src, dst, symlinks=False, dirs_exist_ok=True)
            seen.add(t)
PY

# Collect and compile glib schemas from all bundled packages
SCHEMAS_DIR="${vendor_dir}/glib/share/glib-2.0/schemas"
echo "==> Collecting and compiling GSettings schemas"
# Copy schema files from GTK4 to central location
if [[ -d "${vendor_dir}/gtk4/share/glib-2.0/schemas" ]]; then
  cp "${vendor_dir}/gtk4/share/glib-2.0/schemas"/*.gschema.xml "${SCHEMAS_DIR}/" 2>/dev/null || true
fi
# Compile all schemas
if [[ -d "${SCHEMAS_DIR}" ]] && ls "${SCHEMAS_DIR}"/*.gschema.xml >/dev/null 2>&1; then
  "${BREW_PREFIX}/opt/glib/bin/glib-compile-schemas" "${SCHEMAS_DIR}" 2>&1 | grep -v "does not extend" || true
  echo "  - compiled $(ls ${SCHEMAS_DIR}/*.gschema.xml | wc -l | tr -d ' ') schema files"
else
  echo "  - no schema files found, skipping compilation"
fi

# Generate gdk-pixbuf loaders.cache in bundle
GP_DIR="${vendor_dir}/gdk-pixbuf/lib/gdk-pixbuf-2.0"
if [[ -d "${GP_DIR}" ]]; then
  echo "==> Generating gdk-pixbuf loaders.cache"
  export DYLD_LIBRARY_PATH="${vendor_dir}/gdk-pixbuf/lib:${vendor_dir}/glib/lib:${vendor_dir}/gtk4/lib:${vendor_dir}/pango/lib:${vendor_dir}/cairo/lib:${vendor_dir}/harfbuzz/lib:${vendor_dir}/fribidi/lib"
  mkdir -p "${GP_DIR}/loaders"
  "${BREW_PREFIX}/opt/gdk-pixbuf/bin/gdk-pixbuf-query-loaders" > "${GP_DIR}/loaders.cache" || true
fi

# App launcher
cat > "${APP_MAC}/seatree-launch" <<'LA'
#!/usr/bin/env bash
set -e
APP_DIR="$(cd "$(dirname "$0")"/.. && pwd)"
APP_RES="$APP_DIR/Resources"

# Root for Seatree code
export APPIMAGE_ROOT="$APP_RES/Seatree"
cd "$APPIMAGE_ROOT"

# Brewed runtimes inside app
VND="$APP_RES/vendor/opt"

# Log to user logs so Finder launches are debuggable
LOG_DIR="$HOME/Library/Logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/Seatree.log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "[Seatree] Launching at $(date)"

# Python runtime
if [[ -x "$VND/python@3.13/bin/python3" ]]; then
  PYEXE="$VND/python@3.13/bin/python3"
  export PYTHONHOME="$VND/python@3.13/Frameworks/Python.framework/Versions/3.13"
else
  PYEXE="python3"
fi
# Python packages: python-site contains gi/cairo and other packages
PYTHONPATH_ADD=("$APP_RES/python-site")
# Add Python's standard library site-packages last
_VPY_SITE="$VND/python@3.13/Frameworks/Python.framework/Versions/3.13/lib/python3.13/site-packages"
[[ -d "$_VPY_SITE" ]] && PYTHONPATH_ADD+=("$_VPY_SITE")
export PYTHONPATH="$(IFS=:; echo "${PYTHONPATH_ADD[*]}")"

# Prevent user site-packages from being loaded
export PYTHONNOUSERSITE=1

# Locale safety
export LC_ALL="en_US.UTF-8"
export LANG="en_US.UTF-8"

# Dynamic libs for GTK/GLib etc. and NetCDF
export DYLD_LIBRARY_PATH="$VND/gtk4/lib:$VND/glib/lib:$VND/pango/lib:$VND/cairo/lib:$VND/harfbuzz/lib:$VND/fribidi/lib:$VND/gdk-pixbuf/lib:$VND/gobject-introspection/lib:$VND/netcdf/lib:$VND/netcdf-fortran/lib:$VND/hdf5/lib:$DYLD_LIBRARY_PATH"
export DYLD_FALLBACK_LIBRARY_PATH="$DYLD_LIBRARY_PATH"

# GI runtime
GI_DIRS=()
[[ -d "$VND/gtk4/lib/girepository-1.0" ]] && GI_DIRS+=("$VND/gtk4/lib/girepository-1.0")
[[ -d "$VND/gobject-introspection/lib/girepository-1.0" ]] && GI_DIRS+=("$VND/gobject-introspection/lib/girepository-1.0")
[[ -d "$VND/glib/lib/girepository-1.0" ]] && GI_DIRS+=("$VND/glib/lib/girepository-1.0")
export GI_TYPELIB_PATH="$(IFS=:; echo "${GI_DIRS[*]}"):$GI_TYPELIB_PATH"

# Resources
export GSETTINGS_SCHEMA_DIR="$VND/glib/share/glib-2.0/schemas"
_GP_VER_DIR="$(ls "$VND/gdk-pixbuf/lib/gdk-pixbuf-2.0" 2>/dev/null | head -n1)"
export GDK_PIXBUF_MODULEDIR="$VND/gdk-pixbuf/lib/gdk-pixbuf-2.0/${_GP_VER_DIR}/loaders"
export GDK_PIXBUF_MODULE_FILE="$VND/gdk-pixbuf/lib/gdk-pixbuf-2.0/${_GP_VER_DIR}/loaders.cache"
export GIO_MODULE_DIR="$VND/glib/lib/gio/modules"
export XDG_DATA_DIRS="$VND/gtk4/share:$VND/glib/share:$XDG_DATA_DIRS"

# PATH for helper tools (gawk, gs, bundled GMT)
export PATH="$VND/gawk/bin:$VND/ghostscript/bin:$VND/gmt/bin:$PATH"

echo "[Debug] Python: $PYEXE"
echo "[Debug] PYTHONPATH: $PYTHONPATH"
_GI_SO="$(find "$APP_RES/python-site/gi" -name "_gi*.so" -type f 2>/dev/null | head -1)"
if [[ -n "$_GI_SO" ]]; then
  echo "[Debug] Found _gi extension: $_GI_SO"
else
  echo "[Debug] ERROR: _gi extension not found in $APP_RES/python-site/gi"
fi

exec "$PYEXE" "$APPIMAGE_ROOT/python3/seatree/gui/SEATREE.py"
LA
chmod +x "${APP_MAC}/seatree-launch"

# Minimal Info.plist
cat > "${APP_CNT}/Info.plist" <<'PLIST'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleName</key><string>Seatree</string>
  <key>CFBundleDisplayName</key><string>Seatree</string>
  <key>CFBundleIdentifier</key><string>dev.seatree.Seatree</string>
  <key>CFBundleVersion</key><string>1.0</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleExecutable</key><string>seatree-launch</string>
  <key>CFBundleIconFile</key><string>Seatree</string>
  <key>LSMinimumSystemVersion</key><string>13.0</string>
  <key>NSHighResolutionCapable</key><true/>
</dict>
</plist>
PLIST

# 5) Create DMG (optional)
if [[ "$CREATE_DMG" -eq 1 ]]; then
  mkdir -p "${OUT_DIR}"
  DMG_PATH="${OUT_DIR}/${DMG_NAME}.dmg"
  echo "==> Creating DMG: ${DMG_PATH}"
  TMP_VOL_DIR="${OUT_DIR}/.dmgroot"
  rm -rf "$TMP_VOL_DIR"
  mkdir -p "$TMP_VOL_DIR"
  cp -R "${APP_DIR}" "$TMP_VOL_DIR/"
  hdiutil create -ov -fs HFS+ -srcfolder "$TMP_VOL_DIR" -volname "${VOL_NAME}" -format UDZO "${DMG_PATH}"
  rm -rf "$TMP_VOL_DIR"
  echo "==> DMG created: ${DMG_PATH}"
  echo "Mount test: hdiutil attach \"${DMG_PATH}\""
else
  echo "==> Skipped DMG creation (--no-dmg). App bundle is at ${APP_DIR}."
fi

echo "Done."
