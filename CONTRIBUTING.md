# Contributing

Thank you for contributing to lnc-seeker-hub.

## Development setup

1. Install Rust and Python.
2. Create a local virtual environment if you plan to run the Bokeh UI.
3. Install runtime Python dependencies from `requirements.txt`.
4. Install build tooling from `requirements-build.txt`.
5. Build the Rust extension with `python -m maturin develop`.

## Before opening a pull request

1. Keep changes focused and scoped to one concern.
2. Update documentation when behavior, setup, or configuration changes.
3. Run the relevant Rust checks locally:
   - `cargo check -p lnc-seeker-collect`
   - `cargo test -p lnc_seeker_lib -- --test-threads=1`
4. Do not commit local datasets, generated CSV files, or machine-specific `config.json` values.

## Pull request expectations

1. Describe the user-visible or developer-visible impact.
2. Note any follow-up work or known limitations.
3. Include tests when practical for new logic or bug fixes.

## Reporting issues

Open a GitHub issue for bugs, regressions, or documentation problems. For security-sensitive issues, follow [SECURITY.md](SECURITY.md).