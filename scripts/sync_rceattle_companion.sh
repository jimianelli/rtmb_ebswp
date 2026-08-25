#!/bin/sh
set -eu

source_dir="${RCEATTLE_COMPANION_DIR:-../ebswp_rceattle/docs}"
site_base="https://noaa-afsc.github.io/ebswp_rceattle/"

for required_file in index.html ebswp.pdf
do
  if [ ! -f "$source_dir/$required_file" ]; then
    echo "Missing Rceattle companion artifact: $source_dir/$required_file" >&2
    exit 1
  fi
done

# Keep the historical RTMB URL available while loading all relative report
# resources from the authoritative standalone Rceattle Pages site.
perl -0pe '
  s{((?:href|src)=")(?!#|https?://|mailto:|data:|javascript:)([^"]+)}
   {$1 . "https://noaa-afsc.github.io/ebswp_rceattle/" . $2}ge
' "$source_dir/index.html" > docs/rceattle_ebswp.html

cp "$source_dir/ebswp.pdf" docs/rceattle_ebswp.pdf
