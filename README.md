# Why are common diseases common?

> "Nothing in biology makes sense except in the light of evolution."
> — Theodosius Dobzhansky, 1973

Human diseases vary enormously in prevalence. Some affect a substantial
fraction of the population, whereas others are exceedingly rare. Epidemiology
has documented this variation in extraordinary detail, yet a basic evolutionary
question remains:

**Why are some diseases common while others are rare?**

## Where this question came from

My background is in evolutionary genetics. Earlier in my career, I worked in a
*Drosophila* laboratory, where I was trained in a tradition that places
enormous value on experimental precision: perturb a gene, control the genetic
background, observe the phenotype, and use reverse genetics to establish
mechanism.

That approach is extraordinarily powerful. But I became increasingly
uncomfortable with a view I sometimes encountered: that human genetic studies
are inherently less rigorous or less informative because humans cannot be
experimentally manipulated with the precision possible in model organisms.

My experience in human rare-disease genomics led me to a different view.

Human biology comes with a different kind of scientific resource: an enormous
natural record accumulated through centuries of medicine and, more recently,
through population-scale genomics. Human phenotypes have been observed across
diseases, organs, ages, environments, populations, and generations at a scale
and breadth that no experimental model can reproduce.

Consider Duchenne muscular dystrophy. Generations of physicians, geneticists,
pathologists, molecular biologists, and families have collectively
characterized its incidence, inheritance, natural history, reproductive
consequences, molecular pathology, and clinical progression. It would be
extraordinarily difficult to reproduce that depth of accumulated knowledge for
a single analogous phenotype in a model organism, let alone for thousands of
phenotypes.

This difference shaped how I came to think about human genetics. Experimental
control is one form of scientific power; observing naturally occurring
variation at population scale is another.

There was also a deeper issue that increasingly interested me as an
evolutionary geneticist. A gene knockout can provide exceptionally clean
evidence that a gene affects fitness. But an engineered complete loss of
function is often an extreme point in the genotype–phenotype map, and it may
poorly represent the standing variation responsible for quantitative fitness
differences. The perturbation is real, but it need not represent the
distribution of genetic effects on which evolution ordinarily acts in natural
populations. In that sense, a knockout can resemble the large-effect
loss-of-function mutations I worked on myself in rare human genetic disease,
such as the *DMD* splice-donor variant described below:
extraordinarily informative about biological mechanism, but not necessarily
representative of the genetic variation through which fitness normally varies
within a population.

Evolution in natural populations ultimately unfolds through genetic variation
that arises, segregates, and is transmitted across generations. Much of that
variation is quantitative: many naturally segregating alleles of modest
effect, producing a continuous distribution of phenotypes and fitness among real
individuals. This is what increasingly attracted me to human quantitative
genetics.

A human biobank is messy in ways that a controlled experiment is not. But every
point in its phenotypic distribution corresponds to a real person; every common
allele has actually segregated through human populations; every reproductive
phenotype represents an outcome that actually occurred. The variation is not
experimentally imposed. It is the variation on which human evolution has been
operating.

That gave me a different sense of what counts as evidence in evolutionary
genetics. Three kinds of study ask three different questions:

- **Knockout genetics.** What happens to fitness if I destroy this gene?
- **Rare-disease genetics.** What happens when nature produces an unusually
  severe mutation in this gene?
- **Quantitative population genetics.** How does naturally segregating genetic
  variation relate to fitness across the population?

The third is the one that corresponds to what this project set out to measure.
For the questions I became interested in, that distinction gave me both
confidence and excitement. I was no longer asking only what an experimentally
constructed perturbation could do; I could ask what naturally occurring genetic
variation actually does across a population. I was looking directly at the
distribution of genetic variation, disease, and reproductive success that
exists in a population, with each individual representing one observation from
that evolutionary process.

And behind every data point was a person who once breathed, felt, lived,
reproduced, suffered, and experienced the world much as we do. I found that
realization profoundly moving. These datasets are not merely collections of
measurements; they are the accumulated biological histories of real human
lives. I am deeply grateful to the participants whose contributions made this
collective record of human variation possible.

My earlier work examined how evolutionarily new genes contribute to human
disease:

> Chen J-H, Landback P, Arsala D, Guzzetta A, Xia S, ... Zhang YE, Cheng J,
> Shen B, Long M. Evolutionarily new genes in humans with disease phenotypes
> reveal functional enrichment patterns shaped by adaptive innovation and
> sexual selection. *Genome Research* 35:379–392 (2025).
> doi:10.1101/gr.279498.124

In parallel, I worked directly on rare Mendelian disease, including the
identification of a splice-donor mutation in *DMD* causing Duchenne muscular
dystrophy:

> Chen J, Jia Y, Zhong J, Zhang K, Dai H, He G, Li F, Zeng L, Fan C, Xu H.
> Novel mutation leading to splice donor loss in a conserved site of DMD gene
> causes Duchenne muscular dystrophy with cryptorchidism.
> *Journal of Medical Genetics* 61:741–749 (2024). doi:10.1136/jmg-2024-109896

Rare-disease studies routinely begin with a number: the prevalence of the
disease. Duchenne affects roughly one in several thousand male births; other
disorders may affect one in tens or hundreds of thousands.

Manyuan Long once told me that, in science, finding the right question can be
more important than finding the answer. I kept coming back to that idea. My
intuition was that disease prevalence itself was an important evolutionary
question, not simply a number to report at the beginning of a paper.

That led me to ask something that seemed almost too simple:

What does prevalence mean in evolutionary terms? Why does one genetic disease
remain extremely rare while another complex disease affects millions of people?

## From disease genetics to population genetics

When I began working with Lin Chen and became more deeply involved in
statistical genetics and genomic epidemiology, I returned to this question.

It became apparent that the question could not be answered by experimental
manipulation of a model organism. What made it tractable was precisely the
feature of human genetics that is sometimes regarded as its weakness: we cannot
design the experiment.

Instead, millions of naturally occurring genetic perturbations have already
happened. Recombination has reshuffled them across generations. Natural
selection has acted on them. Physicians have recorded their phenotypic
consequences. Biobanks have measured disease and reproduction in hundreds of
thousands of individuals. GWAS and molecular QTL studies have connected those
phenotypes back to genetic variation.

The experiment, in a sense, has been running for generations. The challenge is
to extract its signal.

This led to a simple evolutionary argument. Natural selection does not act on
disease prevalence itself. Ultimately, it acts through fitness. A
disease-associated allele does not disappear merely because we call it a "risk
allele." Its evolutionary trajectory depends on all of its phenotypic
consequences, including its effects on reproduction.

If genetic effects that increase disease liability are coupled to greater
reproductive success, selection through reproduction may oppose selection
against disease. This is antagonistic pleiotropy. If disease-increasing effects
are instead coupled to reduced reproductive success, the two selective
consequences act in the same direction.

That gives a testable prediction:

**Variation in disease prevalence may partly reflect how disease-associated
genetic variation is coupled to reproductive fitness.**

## Testing the idea

The first stage of this work characterized disease–fertility genetic coupling
across 97 complex diseases. The broader question became testable when
large-scale disease prevalence data could be integrated with genomic estimates
of reproductive effects.

Across thousands of diseases, a systematic pattern emerged. Disease prevalence
increased with more antagonistic disease–fertility coupling: genetic effects
associated with increased disease liability were also associated with greater
reproductive success. Conversely, prevalence decreased with synergistic
coupling, in which disease-increasing effects were associated with reduced
reproductive success.

The same evolutionary logic also led us to examine differences between male and
female reproductive effects, revealing widespread sexually antagonistic effects
among disease-associated genes.

What began as a question raised by rare-disease genetics therefore became a
population-level evolutionary hypothesis:

**The reproductive consequences of disease-associated genetic variation may
help explain why some human diseases are common while others remain rare.**

Long live human evolutionary medicine.

---

Repository layout and per-script notes: [CONTENTS.md](CONTENTS.md)
