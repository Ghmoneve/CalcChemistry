import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Topology.MetricSpace.Basic

open scoped RealInnerProductSpace

example : Continuous fun p : EuclideanSpace ℝ (Fin 3) × EuclideanSpace ℝ (Fin 3) => ⟪p.1, p.2⟫ := continuous_inner
#check @real_inner_comm
#check @real_inner_add_left
#check @real_inner_add_right
#check @isCompact_sphere
#check @ContinuousOn.exists_isMinOn
#check @intermediate_value_Icc
#check @Real.sqrt_le_sqrt
#check @Real.sq_sqrt
#check @ContinuousLinearMap.continuous
