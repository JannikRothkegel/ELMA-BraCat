cmake_minimum_required (VERSION 2.8)

# Standard LeMonADE
add_subdirectory(SimpleSimulatorStandard)

# Branched Catalysis
add_subdirectory(CreatorChainWalkingTertiaryBondWalking)

# Rouse Matrix evaluation of BraCat structures
add_subdirectory(CreatorChainWalkingTertiaryBondWalking_RouseMatrix)

# Topological properties of BraCat structures
add_subdirectory(CreatorChainWalkingTertiaryBondWalking_TopologicalProperties)

# Creational properties of BraCat structures
add_subdirectory(CreatorChainWalkingTertiaryBondWalking_CreationProperties)

# BGRW_CreationProperties
add_subdirectory(BGRW_CreationProperties)

# Rouse Matrix evaluation of BGRW
add_subdirectory(BGRW_RouseMatrix)

# Rouse Matrix evaluation of Hyperstar
add_subdirectory(Hyperstar_RouseMatrix)

# Rouse Matrix evaluation of Hyperdendrimer
add_subdirectory(Hyperdendrimer_RouseMatrix)

