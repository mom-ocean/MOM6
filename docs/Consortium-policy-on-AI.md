# Policy for AI-generated code in MOM6

## Guiding Principle

AI coding assistants are welcome in MOM6 development the same way any other tool is welcome — with the expectation that a human contributor is responsible for everything they submit. AI use is permitted when it supports a contributor's comprehension of the code; it is not permitted as a substitute for that comprehension.


## Contributor Responsibility

Every contribution to MOM6 has a human author who is accountable for its correctness, physical consistency, and conformance to MOM6 coding standards. Every contribution must have a clearly articulated purpose that is succinctly documented in the commit messages and pull request descriptions.
Reviewers may ask any contributor, "Do you understand why this works?" and “Why are you suggesting this?” These questions are always appropriate, and a satisfactory answer is always expected.
Code that the contributor cannot explain should not be submitted.


## Disclosure

Contributors must disclose AI assistance in any pull request or commit where AI tooling materially shaped the submitted code. The tool is not credited as an author; the disclosure is for the benefit of reviewers.
A brief note in the PR description is sufficient, e.g. "Note that an AI tool was used in making this contribution"


## Risks of AI-Assisted Code

The fact that AI-generated code compiles and passes existing tests is not sufficient evidence of correctness. Sign errors, unit mismatches, and flux miscalculations can pass a compiler and look plausible. AI models are capable of producing syntactically correct Fortran; they do not understand ocean physics.

Contributors and reviewers should be alert to the following categories of risk:

Physical correctness: AI tools have no understanding of the underlying science. Parameterizations, numerical schemes, and flux calculations must be verified by a human with domain knowledge, not inferred to be correct because the code runs.
Comprehension debt: If a contributor does not fully understand code they have submitted, that gap compounds over time as the code is further modified, potentially again with AI assistance. See also Section on responsibility.
Vacuous tests: AI-generated tests may pass without actually testing anything meaningful. Tests must be evaluated for genuine coverage and practical relevance.
Code churn: Modifying code without a compelling purpose impedes the process of developing useful capabilities. 
Malicious or unsafe code: AI tools can be manipulated through prompt injection or may reproduce unsafe patterns from their training data. Contributors must review AI-suggested code for suspicious logic, unexpected network calls, file system access, or other behavior inconsistent with the intended change. Any code that cannot be straightforwardly explained by its stated purpose should be treated as suspect.


## Roles for AI in the Workflow

The following uses are encouraged, with appropriate human oversight:

Development assistance: exploring code, planning changes, drafting implementations, explaining unfamiliar code, suggesting refactors.
Debugging: identifying likely error sources, interpreting compiler or runtime messages.
Documentation: drafting docstrings and comments, prose descriptions.
CI triage: understanding feedback or messages on CI failures (e.g., summarizing logs, suggesting fixes for straightforward failures).
Reviewing Pull Requests: AI tools may be used to highlight possible issues in contributions, such as bugs or points of inconsistency with MOM6 practice and policy or clarifying that a contribution is pure refactoring.  However, only a human reviewer may approve a pull request after fully understanding what the contributed code does.

The following uses require particular caution:

Autonomous code contribution with minimal human oversight is not currently accepted.
AI-generated tests must be verified to test something meaningful; vacuously passing trivial tests is a known failure mode.


## Equity

Policy and tooling must not implicitly require paid LLM access. Workflows, templates, and guidance documents must remain useful to contributors without access to commercial AI tools.


## Legal and Institutional Considerations

The IP status of AI-generated code remains unsettled in many countries including the United States, Australia, Korea, the European Union, and the United Kingdom.
Contributors should be aware of their institutional policies before using AI tools.
Contributors are responsible for ensuring that AI tool use complies with their employer's policies.
This consortium does not endorse any specific AI vendor or tool.

