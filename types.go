// Copyright 2012 - 2014 The Seriation Authors. All rights reserved. See the LICENSE file.

package ser

type ObjFn func(sim Matrix64, p IntVector) float64
type Sampler func(mtx Matrix64, sEffort float64) IntMatrix
type OptMethod func(Matrix64, ObjFn, bool) (float64, IntVector)
type OptMethod2 func(Matrix64, IntVector, ObjFn, bool) (float64, IntVector)
type OptMethod3 func(Matrix64, IntVector, ObjFn, bool) float64

type ObjFnQ func(dat Matrix64, pr, pc IntVector) float64
type OptMethodQ func(Matrix64, ObjFnQ, bool) (float64, IntVector, IntVector)
