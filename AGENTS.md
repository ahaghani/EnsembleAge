# AGENTS.md — EnsembleAge

Instructions for any coding agent working in this repository. Every host reads a file by this name
by walking up from the working directory, so this is the single portable place project rules live.

## What this is

An R package for epigenetic age prediction using ensemble clock methods. Public, and intended to be
citable and installable:

```r
devtools::install_github("ahaghani/EnsembleAge")
```

**The package name is load-bearing.** Do not rename the repository or the package — published
citations and cached installs do not follow a redirect.

## ⛔ This repository is PUBLIC

Treat every commit as permanent and world-readable.

- **Nothing internal enters it** — no private paths, no email addresses, no agent tooling
  directories, no editor or assistant configuration, no data, no unpublished results.
- **Nothing internal enters a commit message either**, including co-author or session trailers.
- A path that resolves on one machine is a leak, not a convenience: an absolute path here carries a
  username and often an address. Keep machine-specific settings out of version control entirely.

If tooling needs configuration to work here, configure it **outside this repository**. You cannot
ignore a path without naming it, so the protection is absence rather than an ignore rule.

## Conventions

- Standard R package layout: `R/`, `man/`, `data/`, `inst/`, `vignettes/`, `NAMESPACE`, `DESCRIPTION`.
- `man/` is generated from roxygen comments — edit the comments, not the `.Rd` files.
- Update `NEWS.md` for user-visible changes.
- Keep `R CMD check` clean; the version scheme signals an intent to submit, so a new NOTE is a
  regression.
- Large reference data belongs outside the package, resolved at runtime rather than committed.

## Data

Never write an absolute data path into a script. Where a shared resolver is available, address data
by **logical name** so the code survives the data moving between machines and storage backends.
