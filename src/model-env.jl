# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


# mutable struct Context
#     ndim     ::Int       # Analysis dimension
#     transient::Bool      # Time dependent analysis
#     t        ::Float64   # Time in time dependent analysis
#     outdir   ::String    # Output directory
#     outkey   ::String    # Key name for output files
#     context  ::Context
#     pool     ::ArrayPool
    
#     log      ::IOStream # solver log file
#     alerts   ::IOBuffer # alerts
#     info     ::IOBuffer

#     T      ::Float64  # Pseudo time for current stage
#     ΔT     ::Float64  # Pseudo time current increment
#     Tupdate::Float64  # Pseudo time for last update
#     flushtime::Float64 # Time of last flush
#     residue::Float64  # Ambient absolute temperature in Celsius
#     stage  ::Int      # Current stage
#     inc    ::Int      # Current increment
#     out    ::Int      # Current output file number

#     function Context()
#         this           = new()
        
#         # Analysis related
#         this.ndim      = 3
#         this.transient = false
#         this.t         = 0.0
#         this.outdir    = "."
#         this.outkey    = "out"
#         this.pool      = ArrayPool()

#         this.alerts    = IOBuffer()
#         this.info      = IOBuffer()
        
#         # Stage related:
#         # this.stagebits = StageBits(0.0, 0.0, 0.0, 0, 0, 0, true)
#         this.T       = 0.0
#         this.ΔT      = 0.0
#         this.Tupdate = 0.0
#         this.flushtime = 0.0
#         this.residue = 0.0
#         this.stage   = 0
#         this.inc     = 0
#         this.out     = 0

#         return this
#     end
# end
