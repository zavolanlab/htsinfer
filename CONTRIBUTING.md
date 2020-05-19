# Guidelines for contributing

## General workflow

We are using [Git][res-git], [GitHub][res-github] and [Git Flow][res-git-flow].

> **Note:** If you are a **beginner** and do not have a lot of experience with
> this sort of workflow, please do not feel overwhelmed. We will guide you
> through the process until you feel comfortable using it. And do not worry
> about mistakes either - everybody does them. Often! Our project layout makes
> it very very hard for anyone to cause irreversible harm, so relax, try things
> out, take your time and enjoy the work! :)

## Issue tracker

Please use each project's GitHub [issue tracker][res-issue-tracker] to:

- find issues to work on
- report bugs
- propose features
- discuss future directions

## Submitting issues

Please choose a template when submitting an issue: choose the [Bug report
template][res-bug-report] only when reporting bugs; for all other issues,
choose the [Feature request template][res-bug-report]. Please follow the
instructions in the template.

You do not need to worry about adding labels or milestones for an issue, the
project maintainers will do that for you. However, it is important that all
issues are written concisely, yet with enough detail and with proper
references (links, screenshots, etc.) to allow other contributors to start
working on them. For bug reports, it is essential that they include
reproducible examples.

Please **do not** use the issue tracker to ask usage questions, installation
problems etc., unless they appear to be bugs. For these issues, please use
the [communication channels](#communication) outlined below.

## Communication

If you want to reach out to us (e.g., to discuss potential issues for you to
work on), check the [main documentation][res-docs] for contact information.

## Code style & testing

To make it easier for everyone to maintain, read and contribute to the code,
as well as to ensure that the code base is robust and of high quality, we
would kindly ask you to stick to the following guidelines for code style and
testing.

- Please use a recent version of [Python 3][res-py] (3.7.4+)
- Please try to conform to the used code, docstring and commenting style within
  a project to maintain consistency
- Please use [type hints][py-typing] for all function/method signatures
  (exception: tests)
- Please use the following linters (use default settings unless otherwise
  stated):
  - [`flake8`][res-py-flake8]
  - [`pylint`][res-py-pylint] (use available [configuration][res-py-pylint-conf])
  - [`mypy`][res-py-mypy] OR [`pyright`][res-py-pyright] to help with type hints
- Please use the following test suites:
  - [`pytest`][res-py-pytest]
  - [`coverage`][res-py-coverage]

## Commit messages

In an effort to increase consistency, simplify maintenance and enable automated
change logs, we would like to kindly ask you to write _semantic commit
messages_, as described in the [Conventional Commits
specification][res-conv-commits].

The general structure of _Conventional Commits_ is as follows:

```console
<type>[optional scope]: <description>

[optional body]

[optional footer]
```

Depending on the changes, please use one of the following **type** prefixes:

| Type | Description |
| --- | --- |
| build | The build type (formerly known as chore) is used to identify development changes related to the build system (involving scripts, configurations or tools) and package dependencies.  |
| ci | The ci type is used to identify development changes related to the continuous integration and deployment system - involving scripts, configurations or tools. |
| docs | The docs type is used to identify documentation changes related to the project - whether intended externally for the end users (in case of a library) or internally for the developers. |
| feat | The feat type is used to identify production changes related to new backward-compatible abilities or functionality. |
| fix | The fix type is used to identify production changes related to backward-compatible bug fixes. |
| perf | The perf type is used to identify production changes related to backward-compatible performance improvements. |
| refactor | The refactor type is used to identify development changes related to modifying the codebase, which neither adds a feature nor fixes a bug - such as removing redundant code, simplifying the code, renaming variables, etc. |
| revert | For commits that revert one or more previous commits. |
| style | The style type is used to identify development changes related to styling the codebase, regardless of the meaning - such as indentations, semi-colons, quotes, trailing commas and so on. |
| test | The test type is used to identify development changes related to tests - such as refactoring existing tests or adding new tests. |

In order to ensure that the format of your commit messages adheres to the
Conventional Commits specification and the defined type vocabulary, you can
use the [dedicated linter][res-conv-commits-lint]. More information about
_Conventional Commits_ can also be found in this [blog
post][res-conv-commits-blog].

## Merging your code

Here is a check list that you can follow to make sure that code merges
happen smoothly:

1. [Open an issue](#submitting-issues) _first_ to give other contributors a
   chance to discuss the proposed changes (alternatively: assign yourself
   to one of the existing issues)
2. Clone the repository, create a feature branch off of the default branch
   (never commit changes to protected branches directly) and implement your
   code changes
3. If applicable, update relevant sections of the [documentation][res-docs]
4. Add or update tests; untested code will not be merged; refer to the
   [guidelines](#code-style-&-testing) above for details
5. Ensure that your coding style is in line with the
   [guidelines](#code-style-&-testing) described above
6. Ensure that all tests and linter checks configured in the [Travis
   CI][travis-docs] [continuous integration][ci-cd] (CI) pipeline pass without
   issues
7. If necessary, clean up excessive commits with `git rebase`; cherry-pick and
   merge commits as you see fit; use concise and descriptive commit messages
8. Push your clean, tested and documented feature branch to the remote; make
   sure the [Travis CI][travis-docs] [CI][ci-cd] pipeline passes
9. Issue a pull request against the default branch; follow the instructions in
   the [template][res-pull-request]; importantly, describe your changes in
   detail, yet with concise language, and do not forget to indicate which
   issue(s) the code changes resolve or refer to; assign a project maintainer
   to review your changes

[res-bug-report]: .github/ISSUE_TEMPLATE/bug_report.mdrequest.md
[res-ci-cd]: <https://en.wikipedia.org/wiki/Continuous_integration>
[res-conv-commits]: <https://www.conventionalcommits.org/en/v1.0.0-beta.2/#specification>
[res-conv-commits-blog]: <https://nitayneeman.com/posts/understanding-semantic-commit-messages-using-git-and-angular/>
[res-conv-commits-lint]: <https://github.com/conventional-changelog/commitlint>
[res-docs]: README.md
[res-feature-request]: .github/ISSUE_TEMPLATE/feature_request.md
[res-git]: <https://git-scm.com/>
[res-git-flow]: <https://nvie.com/posts/a-successful-git-branching-model/>
[res-github]: <https://github.com>
[res-issue-tracker]: <https://github.com/zavolanlab/htsinfer/issues>
[res-pull-request]: PULL_REQUEST_TEMPLATE.md
[res-py]: <https://www.python.org/>
[res-py-flake8]: <https://gitlab.com/pycqa/flake8>
[res-py-mypy]: <http://mypy-lang.org/>
[res-py-pylint]: <https://www.pylint.org/>
[res-py-pylint-conf]: pylint.cfg
[res-py-pyright]: <https://github.com/microsoft/pyright>
[res-py-pytest]: <https://docs.pytest.org/en/latest/>
[res-py-coverage]: <https://pypi.org/project/coverage/>
[res-py-typing]: <https://docs.python.org/3/library/typing.html>
[res-travis-docs]: <https://docs.travis-ci.com/>
