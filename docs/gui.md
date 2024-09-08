## GUI Structure

### Project Browser (Left)
- Geometry
  - Multiple Geometry
    - Table (Type, R,theta, Length, Lambda)
    - *Drawing* (2D Sketch)
- Ambient
  - Types (ISA_NOM, ISA_UB, ISA_LB)
  - Type: Custom
    - Table (Altitude vs P,T,Rho)
    - Read From File
- Fluid
    - Types (Perfectgas, hansen)
    - PefectGas
      - Molar Mass
    - Hansen
- Trajectory
  - Type (Velocity, Mach, CFD_DATA, WindTunnel)
  - Velocity
    - Ambient
    - Table (t,Alt,V)
  - Mach
    - Ambient
    - Table (t,Alt,M)
  - CFDDATA
    - Table (t,M,CFD_data)
    - CFD_Data: P,T,M,G,Rho
  - WTunnel
    - Table (t,M,WindTunnelData)
    - WTData: (M,P0,T0)
- BC
  - Radiative BC
    - Radiative Tspace
    - Front/Back
  - Convective BC
    - Front/Back
    - Convection
      - Table (t,h,Tg)
    - Natural convection
      - Table (t,L,Tg)
  - Adiabatic BC
    - (By default)
  - Flux
    - Front/Back
    - Table (t, flux) : Scale factor / Offset
  - Propellant Mass
    - Table (t,PropMass) => mcp/a capacity/Area
    - mass, cp , area (scalar/table)
  - Heat Generation BC
    - Table (t,qdot)

- Material (Local)
  - Option: Inherit from library => options
  - Name
  - Density (scalar/Table) => virgin/cahr/pyrogas
  - Conductivity (scalar/Table) => virgin/cahr/pyrogas
  - Cp (scalar/Table) => virgin/cahr/pyrogas
  - emissivity
  - Habl
  - Tabl
  - Hpyro
  - Tpyro
  - ISUBLM
  - Airnode? option - El : None

- Coordinate System
  - Type: Cartesian, Cylindrical, (Spherical)

- LayerStack (stacking)
  - Collection of Layers
    - Thickness
    - Number of Nodes
    - Material
      - Library/Local

- Solution (Case)
  - Type: Flux
  - Tinitial, Tfinal, tstep (Scalar/Table**)
  - Toutput: Sampling No, Timestep
  - Geometry
  - Trajectory
    - Trajectory (actual: data)
    - h: scale / bias
    - FlowType
    - CpData
    - Shock: On/Off
  - Fluid
  - Params
    - ReC
    - ReTl
  - Options
    - Ttol
    - 

- Solution (Case):
  - Type: Response (2D, Axisym..)
  - CoordSys
  - Tinitial, Tfinal, tstep (Scalar/Table**)
  - Toutput: Sampling No, Timestep
  - Geometry
    - Geometry Data
    - Solution Location: Location (scalar/table) / No of Points => Layerstack
  - Trajectory
    - Trajectory (actual: data)
    - h: scale / bias
    - FlowType
    - CpData
    - Shock: On/Off
  - Fluid
  - BCs
    - Can Add various types
  - __LayerStack ???__ (Temperature dist: overrides IC)
  - ICs
    - Temperature (scalar/Table: Location) => layerstack single value
  - Params
    - ReC
    - ReTl
  - Options
    - Ttol
    - PCtol
    - FlatPlate: Eckert/VD
    - Transition => WettedLength/ Momentum thickness
    - Temperature: units (K/...)

<!-- - Solution (Case):
  - Type: Otimization (2D, Axisym..)
  - CoordSys
  - Tinitial, Tfinal, tstep (Scalar/Table**)
  - Toutput: Sampling No, Timestep
  - Geometry
  - Trajectory
    - Trajectory (actual: data)
    - h: scale / bias
    - FlowType
    - CpData
    - Shock: On/Off
  - Fluid
  - BCs
    - Can Add various types
  - __LayerStack ???__ (Temperature dist: overrides IC)
  - ICs
    - Temperature (scalar/Table: Location) => layerstack single value
  - Params
    - ReC
    - ReTl
  - Options
    - Ttol
    - PCtol
    - FlatPlate: Eckert/VD
    - Transition => WettedLength/ Momentum thickness
    - Temperature: units (K/...)  -->

### Viewports (Centre)
- Geometry
- Results

### Properties (Right)

### Command window (Bottom)