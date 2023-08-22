module PeriodicTableExt

using PeriodicTable: Element, elements

import CrystallographyBase: atomicmass

atomicmass(element::Element) = element.atomic_mass
atomicmass(i::Union{AbstractString,Integer,Symbol}) = elements[i].atomic_mass

end
