# Claude development instructions

These instructions apply to every future coding task in this repository. GitHub is the system of record for completed work. Do not treat a local-only change, an unpushed commit, or a chat response as finalized work.

## Required workflow

1. Before editing, inspect the repository, current branch, working-tree status, existing documentation, and relevant tests. Confirm the requested scope and acceptance criteria. Do not overwrite or discard existing work.
2. Create a focused task branch from the appropriate up-to-date base branch. Never force-push or commit directly to the protected/shared main branch unless the user explicitly authorizes an exception.
3. Implement the requested change, including appropriate tests and updates to documentation. Keep commits focused and descriptive.
4. Run relevant tests, formatting, linting, type checks, and data-validation checks where available. Report exactly what ran and any failures. Never claim unexecuted tests passed.
5. Review the diff for unintended changes, secrets, credentials, personal or sensitive data, and large generated artifacts. Do not commit these. Use .gitignore and safe sample data where appropriate.
6. Commit the intended files, push the task branch to GitHub, and open or update a pull request targeting the correct base branch. The PR should describe purpose, changes, testing, risks, and remaining work.
7. Verify that the pushed commit and pull request exist on GitHub. Report the branch name, commit SHA, PR URL, test results, and whether review/merge is pending.
8. Only mark a task as finalized when its work is durably pushed to GitHub and its PR is available. Mark it **pending review/merge** until merged. If pushing or creating a PR is blocked by authentication, permissions, connectivity, or missing tools, state **NOT FINALIZED**, preserve the work, and provide the exact next steps; never imply that GitHub was updated.
9. Do not merge a PR, delete remote branches, change repository visibility, alter CI/security protections, publish releases, or deploy to production without explicit authorization. Never bypass failing checks to claim completion.

## Collaboration with ChatGPT

- GitHub issues, branches, pull requests, and repository documentation are the shared source of truth across Claude and ChatGPT.
- For substantial work, record the goal and acceptance criteria in the issue or PR so another assistant can review without relying on chat history.
- Make PRs small enough for independent review; provide a concise handoff listing changed files, design decisions, commands run, test results, and unresolved questions.
- Incorporate independent review feedback in follow-up commits on the same PR; push and re-run relevant checks.

## Scientific and data-analysis work

- Document inputs, sample metadata, references and versions, parameters, dependencies, random seeds where relevant, and output locations.
- Prefer reproducible scripts and configuration over undocumented manual notebook steps.
- Validate inputs and include small synthetic or appropriately shareable test fixtures. Do not commit restricted datasets, identifying research information, secrets, or bulky derived files.
- Clearly distinguish exploratory findings from validated conclusions and document statistical assumptions and limitations.

## Final response template

**GitHub status:** Pushed / Pending review / Merged / NOT FINALIZED
**Repository and branch:** ...
**Commit SHA:** ...
**Pull request:** ...
**Changes:** ...
**Checks run and results:** ...
**Known limitations or blockers:** ...
**Next action:** ...
