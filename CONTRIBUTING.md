## Contribution guidelines

When you do a pull request, please open it to the `dev` branch (not `master`).

Each time you make a changes in the code, please write a note in the [CHANGELOG.md](CHANGELOG.md) file. 
We follow the format from [KeepChangelog.com](https://keepachangelog.com/en/). So, if it's a change in the `dev` branch, 
please write about your changes in `## [Unreleased]` section. Then, when the changes are merged to `master`, the title is 
changed to the corresponding version number.

If you merge any branch to `master`, please follow the next steps:
1. Run "walkthrough.Rmd" to ensure that the changes didn't break any functionality
2. Increase version of the package in the [DESCRIPTION](DESCRIPTION) file. We use [Semantic Versioning](https://semver.org/) guidelines. Briefly:
    - all small fixes increase the 3'rd digit
    - fixes, which change old behavior, i.e. break reproducibility increase 2'nd digit
    - all large changes or changes, which break the user interface increase the 1'st digit
3. Ensure that [CHANGELOG](CHANGELOG.md) contains information about all changes. Update the upper title from `## [Unreleased]` to the corresponding version number.
4. Tag the version for git release and docker rebuild:
    - `git tag -a v1.1.2 -m "v1.1.2"`
    - `git push origin master --tags`
    - Alternatively, after pushing changes to master, on GitHub open "Releases" - "[Draft new release](https://github.com/hms-dbmi/conos/releases/new)", fill the tag name and description using the latest changes from CHANGELOG
