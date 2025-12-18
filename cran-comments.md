## Resubmission

This is a resubmission, with the following changes made to address maintainer
feedback:

* There are currently no references for the novel methods implemented here (a
  methods paper is in preparation), so this resubmission does not add 
  references to the DESCRIPTION
* Limited all examples to maximum 2 cores, per CRAN policy
* Avoided all use of \dontrun{} in examples
* Avoid modifying .GlobalEnv by changing with_seed() utility

All checks pass cleanly.
