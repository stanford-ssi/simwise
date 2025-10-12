# Testing

For each folder of tests, create a `runtests.jl` file (name really doesn't matter as long as it is included in the parent folder's `runtests.jl`), and include all subfolders and files. Also name the tests appropriately. Verbose is better.

Inside the unit tests, replicate each `src/` folder in this folder so that we can figure out which tests are for which cases.

Plan for later:
- Including higher level `integration_test` folders for scenarios and everything else