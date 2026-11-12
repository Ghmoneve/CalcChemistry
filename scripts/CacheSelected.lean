import Cache.Requests

/-!
Targeted transport fallback for pinned mathlib sources. Uses the upstream
hashing, downloading and unpacking functions unchanged, but does not request
optional ProofWidgets JavaScript via Git tags (archive installs lack .git).
This changes delivery only, never theorem statements or kernel checking.
Run with `lake env lean --run scripts/CacheSelected.lean Module.Name ...`.
-/
open Cache Cache.IO Cache.Hashing Cache.Requests System

def main (args : List String) : IO Unit := do
  CacheM.run do
    checkForToolchainMismatch
    let roots ← parseArgs ("get" :: args)
    if roots.isEmpty then
      throw <| IO.userError "Pass explicit mathlib module names"
    let memo ← getHashMemo roots
    let hashes ← memo.filterByRootModules roots.keys
    validateLeanTar
    downloadFiles hashes false true
    unpackCache hashes false
