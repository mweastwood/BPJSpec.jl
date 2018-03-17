# Enormous Matrices

Fundamentally, $m$-mode analysis involves operations on block-diagonal matrices. Although the
block-diagonal structure of these matrices makes it possible to store and perform computations on
these matrices, they can still be enormous. That is, while each block individually may fit in the
system memory, multiple blocks may not. Here we introduce some of the abstractions developed for
working with block diagonal matrices that may each be multiple terabytes in size.

## API

```@meta
CurrentModule = BPJSpec
DocTestSetup = quote
    using BPJSpec
    using Unitful
end
```

```@docs
AbstractBlockMatrix
SimpleBlockVector
SimpleBlockMatrix
MBlockVector
MBlockMatrix
FBlockVector
FBlockMatrix
MFBlockVector
MFBlockMatrix
LBlockMatrix
LMBlockVector
```

