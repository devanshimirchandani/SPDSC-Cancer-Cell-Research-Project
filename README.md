# Cancer-Cell-SPDSC-Research-Project

## Table of Contents

- [Overview](#overview)
- [Contributing](#contributing)
	- [Commit Messages](#commit-messages)
	- [Branch Names](#branch-names)
	- [Pull Requests](#pull-requests)
		- [Title](#title)
		- [Description](#description)
		- [Example Pull Request](#example-pull-request)
- [Developing Team](#developing-team)
- [Project Resources](#project-resources)

## Overview

Exploring predictive modelling of disease states based on cell type distributions in spatial transcriptomics data.

## Contributing

Read the guidelines below to write good commit messages, and branch names, and make pull requests that follow the conventions we will be using throughout the project.

### Commit Messages

- Capitalise the subject line.
- Do not end with a period.
- Use imperative mood, i.e. instead of *"Added ..."* write *"Add ..."*.
- Keep messages logical and relevant, do not write things like *"Please work"* or *"I hate frontend"*. To help decide the extent of this, imagine trying to access a point in the project 2 weeks ago, it would be better to have something like *"Add CSS for Navbar template"* or , so that we know from a glance what the commit is for.
- For more detailed messages, use `git commit -m <title> -m <description>`, however short and concise is still preferred.

### Branch Names
Make a branch using `git checkout -b <branch_name>`.
- Names fall under one of **4** categories
	- Minor Feature: `minor-FeatureName`
	- Major Feature: `major-FeatureName`
	- Patch: `patch-PatchName`
	- Miscellaneous: `name`
		- For example `documentation` for changing the README, or adding another markdown

### Pull Requests
*Summarised from [this article](https://namingconvention.org/git/pull-request-naming.html).*

#### Title
- Short and descriptive summary
- Start with corresponding ticket/story id (e.g. from Jira, GitHub issue)
- Should be capitalized and written in imperative present tense
- Do not end with period
- Suggested format: *#<Ticket_ID> PR description*

#### Description
- Separated with a blank line from the subject
- Explain what, why, etc.
- Max 72 chars
- Each paragraph capitalized

#### Example Pull Request
```
This pull request is part of the work to make it easier for people to contribute to naming convention guides. One of the easiest way to make small changes would be using the Edit on Github button.

To achieve this, we needed to:
- Find the best Gitbook plugin which can do the work
- Integrate it in all the pages to redirect the user to the right page on GitHub for editing
- Make it visible on the page so users can notice it easily
```

## Developing Team
- [Devanshi Mirchandani](https://github.com/devanshimirchandani)
- Nancy [Mnasnan Seamorntham](https://github.com/mnasnan)
- [Veronica](https://github.com/Veronicazwt)

## Project Resources
- [Seurat analysis tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
- [Seurat command list](https://satijalab.org/seurat/articles/essential_commands.html)
- [Colon cancer paper](https://www.nature.com/articles/s41588-022-01088-x)
