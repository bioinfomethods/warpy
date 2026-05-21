# Bpipe Pipeline Framewok

Bpipe is a pipeline framework that converts regular Groovy code containing
closure definitions into pipeline stages, implementing an idiomatic
data flow between pipeline stages with built in dependency tracking.

## Core idiom

If a closure is declared as:

```
foo = {
    exec """
        cp -v $input $output
    """
}
```

Then within Bpipe, `foo` automatically becomes defined as a "pipeline stage" that
is subjected to managed execution, for example, it may be sent to a cloud system or an
HPC cluster job.

Within `exec` statements, contents are executed as `bash` scripts, with some important modifications:

- adjacent lines are JOINED together. So there is no need to backslash terminate lines that should be part
  of the same command: that will happen automatically
- special variables are provided: `input`, `output` and `branch`. The input and output variables
  are selectors that map to dynamically built file names. They   
  have special behaviour for file extensions: the file extension acts as a filter on the file
  names that are built. The `input` names map to upstream
  files - they search backwards in the pipeline flow to find the first file that matches the
  file extension. The `output` file names are built by appending the stage name and any suffix
  provided to the input file name. This behaviour is modified by `transform` or `produce`
  wrappers that change the output file names built to replace file extensions or replace the
  whole file name respectively
- Pipeline flow: pipeline stages (defined as closures) are joined using operators:
  - the + operator for closures
    is implemented to execute stages sequentially, passing the outputs of the first stage as the inputs to the next
    Example:
    ```
    // Runs stages sequentially flowing outputs from earlier stages to later stages as inputs
    run {  foo + bar + baz }
    ```
  - lists of stages execute in parallel:
    ```
    // baz+bop and bar run in parallel
    run { foo + [ baz + bop, bar ] }
    ```
  - a list multiplied by a list runs the elements in parallel in the second list for every element of the first list:
    ```
    // Run foo and bar both for 'a','b','c'
    run { ['a','b','c'] * [ foo, bar ] }
    ```
    Inside the parallel branches, the branch has a name assigned from the value in the list (a,b,c above), accessible as
    `branch.name`
- Input resolution errors: referencing an input that does not exist (e.g. `input[filetype]` where no file
  matches the given extension) will cause Bpipe to throw an error, not return null. This means that any
  reference to `input[x]` — even in log statements or println calls — must use a value of `x` that
  resolves to an actual file in the pipeline flow. Be careful to distinguish between variables used for
  input resolution (which must match real file extensions) and variables used for metadata/API purposes
  (which may differ from the actual file extension).
-  Settinv variable values with "using" syntax: Closures have a method added to their metaClass that adds variables to the binding of the closure via
  `using`:
  ```
  // execute baz with the value of the `a` variable defined as 1
  run { foo + baz.using(a:1) 
  ```
  
## Branches and Branch Variables  

- Branch variables: parallel running segments form a branch. A branch variable is defined using `branch.<variable> = <value>`
  - these variables flow to all steps in the branch, including any nested parallel segments (branches)
- File path resolution: Do NOT manually construct file paths (e.g. `"results/" + run_id + '_' + sample + '.summary.karyotype.tsv'`).
  Instead, use Bpipe's input resolution to locate files by their extension. This ensures the pipeline
  remains robust to changes in directory structure or naming conventions. For example:
  ```
  // BAD - fragile, breaks if paths change:
  def karyotype_file = "results/" + run_id + '_' + sample + '.summary.karyotype.tsv'
  def lines = new File(karyotype_file).readLines()

  // GOOD - uses Bpipe input resolution:
  def lines = new File(input.karyotype.tsv.toString()).readLines()
  ```
  Bpipe will automatically resolve the correct file for the current sample/branch context.
  
## Stage Variables

Stage variables are defined using a special custom Bpipe construct using a `var` keyword. The `var` is distinct from Java usage, and defines
optional variables in the form

```groovy
foo = {
    var bar : 'hello'
    ...
}
```

In this example, the variable bar is assigned the default of 'hello' only if it does not already have a value. Prefer to use
this method to define optionally overriden variables within stages.
