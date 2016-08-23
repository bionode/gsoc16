# GSoC Bionode 2016

For my (Julian Mazzitelli, @thejmazz) GSoC project, I built [bionode-watermill](https://github.com/bionode/bionode-watermill) under [bionode](https://bionode.io) as part of the [Open Bioinformatics Foundation](https://www.open-bio.org/wiki/Main_Page). I can't thank enough the support of my mentors [bmpvieira](https://github.com/bmpvieira), [yannickwurm](https://github.com/yannickwurm), [mafintosh](https://github.com/mafintosh) and [maxogden](https://github.com/maxogden).

[![gitter][gitter-image]][gitter-url]

## Workflow engine

As of August 23, 2016, I am the only contributor to the project. List of commits [here](https://github.com/bionode/bionode-watermill/commits?author=thejmazz). The project readme describes the current status and the remaining roadmap in greater detail. 

This project was [originally described](https://summerofcode.withgoogle.com/projects/#6585953724399616) as a "workflow engine for streamed data analysis" with the features:

- organize tasks into streams
- allow parallelization of streams
- work locally as well as with common HPC systems
- integrate [dat](http://dat-data.com/)

The project plan took on a slightly new direction after I had time to properly investigate **the problem in bioinformatics workflow systems**. As a result, **streams became a secondary goal of the project**.

My investigation into the alternative workflow tools (bash, make, snakemake, nextflow) is detailed here in this [blog post](https://github.com/thejmazz/jmazz.me/blob/a7391bd385223dc9b9ee20f16ac0dc7122cefd65/content/post/ngs-workflow.md). I realized that a main concern of these tools is the handling of **input and output files**. While the channels of Nextflow improve the quality of describing a workflow over prerequisites and targets from make/snakemake, there was still room for improvement in my opinion. One of the specific goals of my project was to enable the **reuse of the same task in different places**. For example, if `taskA1` produces a file that matches `*.bam`, and `taskA2` does the same (though by different means), one should be able to use the same `taskB` which looks for input from a `*.bam` as the followup to *both* `taskA1` and `taskA2`. The existing alternative is to produce file patterns that identify where they came from (e.g. `*.a1.bam` and `*.a2.bam`, or two Nextflow channels) and two versions of `taskB`, each to handle only a slightly different input. **This seriously impedes the experimentation process when creating scientific workflows**. The blog post [NGS Workflows](https://github.com/thejmazz/jmazz.me/blob/a7391bd385223dc9b9ee20f16ac0dc7122cefd65/content/post/ngs-workflow.md) ends with a proposed JavaScript implementation. **I have spent the first third of GSoC coming up with that proposal as a solution to existing issues, and the last two thirds implementing an MVP for that same proposal**. With bionode-watermill, given that a set of task declarations exist, one can describe a worfklow like so:

```javascript
const pipeline = join(
  junction(
    join(getReference, bwaIndex),
    join(getSamples, fastqDump)
  ),
  trim, mergeTrimEnds,
  decompressReference, // only b/c mpileup did not like fna.gz
  join(
    fork(filterKMC, filterKHMER),
    alignAndSort, samtoolsIndex, mpileupAndCall // 2 instances each of these
  )
)

// Promise
pipeline().then(function (results) { ... } ).catch(function (err) { ... })
// Or callback
pipeline(function (err, data) { ... })
```

For the sake of context, a task is defined like so:

```javascript
const mpileupAndCall = task({
  input: {
    reference: '*_genomic.fna',
    bam: '*.bam',
    _bai: '*.bam.bai'
  },
  output: 'variants.vcf',
  name: 'samtools mpileup | bcftools call',
}, ({ input }) => `
samtools mpileup -uf ${input.reference} ${input.bam} | \
bcftools call -c - > variants.vcf
`)
```

See the [full example](https://github.com/bionode/bionode-watermill/blob/d2dba9f16066da1438a7f3d5c72381dd249f4b4f/examples/variant-calling-filtered/pipeline.js). Here, each task is completely swappable with another as long as the input and output patterns are the same. The `join` and `junction` functions allow the user to describe the flow of the whole pipeline in one place, given that each task is defined above. This is an improvement over methods where pipeline flow is described via properties within task declarations (e.g. channels) and lets the user witness the whole pipeline in one place. The `fork` function allows the user to describe a set of tasks which **take the same input, and produce the same output, yet through different means, as one atomic piece of the pipeline**. All tasks following the fork will automatically be instanced for each forking task. This could be used for example to:

- compare two different programs
- compare paramaters for the same program

*Effortlessly* - simply wrap the task with `fork`.

## What is Left To Be Done?

Currently, bionode-watermill **is not a drop in replacement for snakemake or nextflow**. It is an **MVP illustrating how to compose pipelines where each task describes input/ouput using file glob patterns**. Furthermore, **streams can only be used internally in tasks, and not between tasks**. This is far from the orginal "streaming pipeline system". However, there is a roadmap for all these features, and the code has been structured in a way that (once you know redux basics) it is simple to follow the flow of control throughout the codebase. Some features that should be implemented to improve the readyness of bionode-watermill for "production use":

- streaming between tasks: The plans to implement this are: if `input` or `output` are not defined on a task, `join` will recognize that and pipe output from one task into input for the other task. `fork` will have to use things like [multi-write-stream](https://github.com/mafintosh/multi-write-stream). Use case: have a search task that queries NCBI for SRA IDs, pipe that into a filtering task, pipe that into a task that triggers a pipeline for each SRA ID.
- metrics: real time, cpu time, memory usage, file sizes, etc. On an individual task basis, for joins, for the whole pipeline. 
- UI, DAG visualization: a complex DAG is produced for the pipeline and logs can be emitted from different tasks at once. A web/electron app that can be used to view this in real time would help to grasp the pipeline. A graph visualization is much easier to understand than a JSON dump. As well, since redux was used, maintaining state between client/server should be *relatively manageable*. Going further, if the UI could edit/create tasks, that would be very cool.
- Example pipelines from published papers. This is important to gain recognition from the bioinformatics community, find use cases and edge cases, develop usage patterns. Current examples can be found at [bionode-watermill/examples](https://github.com/bionode/bionode-watermill/tree/d2dba9f16066da1438a7f3d5c72381dd249f4b4f/examples).
- Improved file validation: custom validators described from task level. Use case: pass a function, that given two paired reads, checks if each has same number of spots. Currently, the file validation checks for non null files, and that the hash of a previously ran task's output matches the current file when checking for resumabiltility (i.e. do not skip a task if the output file has a different hash than what was recorded previously).
- Improve stability of `fork`. Right now it works fine in some cases but may break in others. `fork` is only tested in the case where it is within a join, and followed by some number of tasks (but not `junction` or `join`) - e.g. `join(fork(A1, A2), t1, t2, ... , tn))`
- Modularize the codebase into smaller scoped, well tested modules. The code is modular already but there is still room for improvement. Operations like `join`, `junction`, `fork`, for example could be generalized into their own "orchestration" module and the task lifecycle management could be isolated into a module which exposes an API such as `lifecycle.addStage()`, `stage.onSuccess(anotherStage)` - almost like a state machine. The stability of the lifecycle management, management of context sharing between orchestration operations, passing of context into a new task, composability of orchestration operators are some of the things that need to strictly defined, well tested, and behave well under errors. There is a lot of moving parts in this application, and it can become difficult to keep everything in check. Redux helps as a pattern to make when/how changes to the state occur, but may or may not be as necessary moving forward.
- Much more: see the [project issues](https://github.com/bionode/bionode-waterwheel/issues). I am available to provide feedback on any questions, suggest possible implementation in detail, etc. Any activity in the issues will make me very happy!


[gitter-image]: https://img.shields.io/gitter/room/bionode/gsoc16.svg?style=flat-square
[gitter-url]: https://gitter.im/bionode/gsoc16
