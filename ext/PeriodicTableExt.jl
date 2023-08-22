module PeriodicTableExt

using CrystallographBase
using PeriodicTable: Element, elements

import CrystallographBase: atomicmass

atomicmass(element::Element) = element.atomic_mass
atomicmass(i::Union{AbstractString,Integer,Symbol}) = elements[i].atomic_mass

end
