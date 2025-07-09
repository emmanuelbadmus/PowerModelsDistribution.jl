const _phasecode_map = Dict(
    "PhaseCode.ABCN" => [1, 2, 3],
    "PhaseCode.ABC" => [1, 2, 3],
    "PhaseCode.AB" => [1, 2],
    "PhaseCode.AC" => [1, 3],
    "PhaseCode.BC" => [2, 3],
    "PhaseCode.A" => [1],
    "PhaseCode.B" => [2],
    "PhaseCode.C" => [3],
    "PhaseCode.AN" => [1],
    "PhaseCode.BN" => [2],
    "PhaseCode.CN" => [3]
)

_phase_map = Dict(
    "SinglePhaseKind.A" => 1,
    "SinglePhaseKind.B" => 2,
    "SinglePhaseKind.C" => 3,
    "SinglePhaseKind.N" => 4
)

const _multipliers_map = Dict(
    "UnitMultiplier.m" => 1e-3,
    "UnitMultiplier.c" => 1e-2,
    "UnitMultiplier.d" => 1e-1,
    "UnitMultiplier.da" => 1e1,
    "UnitMultiplier.h" => 1e2,
    "UnitMultiplier.k" => 1e3,
    "UnitMultiplier.M" => 1e6,
    "UnitMultiplier.G" => 1e9,
    "UnitMultiplier.T" => 1e12,
    "UnitMultiplier.P" => 1e15,
    "UnitMultiplier.E" => 1e18,
    "UnitMultiplier.Z" => 1e21
)


"initializes the base math object of any type"
function _init_math_obj_ravens(obj_type::String, eng_id::Any, eng_obj::Dict{String,<:Any}, index::Int; pass_props::Vector{String}=String[])::Dict{String,Any}
    math_obj = Dict{String,Any}(
        "name" => "$eng_id",
        "source_id" => "$obj_type.$eng_id"
    )

    math_obj["index"] = index

    return math_obj
end


"converts impendance in Ohm/m by multiplying by length"
function _impedance_conversion_ravens(eng_obj::Dict{String,<:Any}, vals::Matrix{Float64})
    return vals .* get(eng_obj, "Conductor.length", 1.0)
end


"converts admittance by multiplying by 2πωl"
function _admittance_conversion_ravens(eng_obj::Dict{String,<:Any}, vals::Matrix{Float64})
    2.0 .* pi .* vals .* get(eng_obj, "Conductor.length", 1.0) ./ 1e9
end


"converts impendance in Ohm/m by multiplying by length"
function _impedance_conversion_ravens(data_eng::Dict{String,Any}, eng_obj::Dict{String,Any}, key::String)

    _conductor_count =  data_eng["PerLengthPhaseImpedance.conductorCount"]
    _impedance_matrix = zeros(Float64, _conductor_count, _conductor_count)

    for obj in data_eng["PerLengthPhaseImpedance.PhaseImpedanceData"]
        row = obj["PhaseImpedanceData.row"]
        col = obj["PhaseImpedanceData.column"]
        value = get(obj, key, 0.0)
        _impedance_matrix[row, col] = value
        _impedance_matrix[col, row] = value
    end

    return _impedance_matrix .* get(eng_obj, "Conductor.length", 1.0)
end


"converts impendance in Ohm/m by multiplying by length"
function _impedance_conversion_ravens_energy_source(data_eng::Dict{String,Any}, eng_obj::Dict{String,Any}, key1::String, key2::String)
    # Default energy sources considered 3 phases
    nphases = 3
    _impedance_matrix = zeros(Float64, nphases, nphases)

    z = get(eng_obj, key1, 0.0)
    z0 = get(eng_obj, key2, 0.0)

    for i in 1:nphases
        for j in 1:i
            if(i==j)
                _impedance_matrix[i, j] =  z + ((z0 - z)/3)
            else
                _impedance_matrix[i, j] = (z0 - z)/3
                _impedance_matrix[j, i] = (z0 - z)/3
            end
        end
    end

    return _impedance_matrix .* get(eng_obj, "Conductor.length", 1.0)
end


"converts admittance by multiplying by 2πωl"
function _admittance_conversion_ravens(data_eng::Dict{String,<:Any}, eng_obj::Dict{String,<:Any}, key::String)

    _conductor_count =  data_eng["PerLengthPhaseImpedance.conductorCount"]
    _admittance_matrix = zeros(Float64, _conductor_count, _conductor_count)

    for obj in data_eng["PerLengthPhaseImpedance.PhaseImpedanceData"]
        row = obj["PhaseImpedanceData.row"]
        col = obj["PhaseImpedanceData.column"]
        value = get(obj, key, 0.0)
        _admittance_matrix[row, col] = value
        _admittance_matrix[col, row] = value
    end

    return _admittance_matrix .* get(eng_obj, "Conductor.length", 1.0) ./ 2 # divide by 2 to get both sides _to and _fr
end

"extracts the name from a ravens reference string"
function _extract_name(element)

    name = replace(split(element, "::")[2], "'" => "")
    return name
end


"calculates the shunt admittance matrix based on terminals and b or g"
function _calc_shunt_admittance_matrix(terminals, b)

    N = length(terminals)
    _shunt_matrix = b* Matrix(LinearAlgebra.I, N, N)
    return _shunt_matrix

end


"""
    apply_voltage_bounds_math!(data::Dict{String,<:Any}; vm_lb::Union{Real,Missing}=0.9, vm_ub::Union{Real,Missing}=1.1)

add voltage bounds to all buses based on per-unit upper (`vm_ub`) and lower (`vm_lb`), in MATHEMATICAL.
"""
function apply_voltage_bounds_math!(data::Dict{String,<:Any}; vm_lb::Union{Real,Missing}=0.9, vm_ub::Union{Real,Missing}=1.1)
    if ismultinetwork(data)
        for (_, nw_data) in data["nw"]
            for (_, bus_data) in nw_data["bus"]
                if (bus_data["bus_type"] != 3)
                    bus_data["vmin"] = ones(length(bus_data["vmin"])).*vm_lb
                    bus_data["vmax"] = ones(length(bus_data["vmax"])).*vm_ub
                end
            end
        end
    else
        for (_, bus_data) in data["bus"]
            if (bus_data["bus_type"] != 3)
                bus_data["vmin"] = ones(length(bus_data["vmin"])).*vm_lb
                bus_data["vmax"] = ones(length(bus_data["vmax"])).*vm_ub
            end
        end
    end
end



function build_base_voltage_graphs(data::Dict{String,<:Any})::Tuple{Dict{Int,String}, Graphs.SimpleGraph}
    nodes = Dict(cn => n for (n, cn) in enumerate(keys(data["ConnectivityNode"])))
    G = Graphs.SimpleGraph(length(nodes))

    for edge_dicts in [
        _recursive_dict_get(data, ["PowerSystemResource", "Equipment", "ConductingEquipment", "Conductor", "ACLineSegment"], Dict()),
        _recursive_dict_get(data, ["PowerSystemResource", "Equipment", "ConductingEquipment", "Switch"], Dict())
    ]
        for (i, edge) in edge_dicts
            Graphs.add_edge!(G, nodes[match(Regex("ConnectivityNode::'(.+)'"), edge["ConductingEquipment.Terminals"][1]["Terminal.ConnectivityNode"]).captures[1]], nodes[match(Regex("ConnectivityNode::'(.+)'"), edge["ConductingEquipment.Terminals"][2]["Terminal.ConnectivityNode"]).captures[1]])
        end
    end

    return Dict{Int,String}(n => cn for (cn, n) in nodes), G
end

function find_voltages(data::Dict{String,<:Any})::Dict{String,Any}
    voltages = Dict{String,Any}()

    for (i, es) in _recursive_dict_get(data, ["PowerSystemResource", "Equipment", "ConductingEquipment", "EnergyConnection", "EnergySource"], Dict())
        nominal_voltage = get(es, "EnergySource.nominalVoltage", missing)
        if !ismissing(nominal_voltage)
            voltages[match(Regex("ConnectivityNode::'(.+)'"), es["ConductingEquipment.Terminals"][1]["Terminal.ConnectivityNode"]).captures[1]] = nominal_voltage
        end
    end

    for (i, tr) in _recursive_dict_get(data, ["PowerSystemResource", "Equipment", "ConductingEquipment", "PowerTransformer"], Dict())
        info_name = match(Regex("TransformerTankInfo::'(.*)'"), get(get(tr, "PowerTransformer.TransformerTank", [Dict()])[1], "PowerSystemResource.AssetDatasheet", "TransformerTankInfo::''")).captures[1]
        trinfos = _recursive_dict_get(data, ["AssetInfo", "PowerTransformerInfo", info_name, "PowerTransformerInfo.TransformerTankInfos", info_name, "TransformerTankInfo.TransformerEndInfos"], [])
rated_u = merge(
            filter(x -> !ismissing(x.second), Dict(get(trinfo, "TransformerEndInfo.endNumber", n) => get(trinfo, "TransformerEndInfo.ratedU", missing) for (n, trinfo) in enumerate(trinfos))),
            filter(x -> !ismissing(x.second), Dict(get(pte, "endNumber", n) => get(pte, "PowerTransformerEnd.ratedU", missing) for (n, pte) in enumerate(get(tr, "PowerTransformer.PowerTransformerEnd", []))))
        )
        for (n, term) in enumerate(get(tr, "ConductingEquipment.Terminals", []))
            voltages[match(Regex("ConnectivityNode::'(.+)'"), term["Terminal.ConnectivityNode"]).captures[1]] = get(rated_u, n, missing)
        end
    end

    return voltages
end


function find_base_voltages(data::Dict{String,<:Any})::Dict{String,Any}
    node_lookup, G = build_base_voltage_graphs(data)

    voltages = find_voltages(data)

    ccs = Graphs.connected_components(G)

    voltage_per_cc = Dict()
    for (n, cc) in enumerate(ccs)
        for i in cc
            if node_lookup[i] in keys(voltages)
                voltage_per_cc[n] = voltages[node_lookup[i]]
                break
            end
        end
    end

    return Dict{String,Any}(node_lookup[i] => get(voltage_per_cc, n, missing) for (n, cc) in enumerate(ccs) for i in cc)
end

function _recursive_dict_get(dict::Dict, path::Vector{<:Any}, default::Any)::Any
    if length(path) > 1
        return _recursive_dict_get(get(dict, path[1], Dict()), path[2:end], default)
    else
        return get(dict, path[1], default)
    end
end

function _recursive_dict_set!(dict::Dict, path::Vector{<:Any}, value::Any)
    if length(path) > 1
        _recursive_dict_set!(dict[path[1]], path[2:end], value)
    else
        dict[path[1]] = value
    end
end

function add_base_voltages!(data::Dict{String,<:Any}; overwrite::Bool=false)::Nothing
    if overwrite || "BaseVoltage" ∉ keys(data)
        data["BaseVoltage"] = Dict{String,Any}()
    end

    base_voltages = find_base_voltages(data)

    unique_bv = unique(values(base_voltages))

    for bv in unique_bv
        data["BaseVoltage"]["RAVENS_BaseV_$(bv/1000.0)_kV"] = Dict{String,Any}(
            "IdentifedObject.name" => "RAVENS_BaseV_$(bv) V",
            "IdentifedObject.mRID" => "$(UUIDs.uuid4())",
            "Ravens.cimObjectType" => "BaseVoltage",
            "BaseVoltage.nominalVoltage" => bv
        )
    end

    encon_path = ["PowerSystemResource", "Equipment", "ConductingEquipment", "EnergyConnection"]

    for path in [["EnergyConsumer"], ["EnergySource"], ["RegulatingCondEq", "PowerElectronicsConnection"], ["RegulatingCondEq", "RotatingMachine"]]
        path = [encon_path..., path...]
        for (i, item) in _recursive_dict_get(data, path, Dict())
            if !overwrite && "ConductingEquipment.BaseVoltage" ∈ keys(item)
                continue
            else
                cn = match(Regex("ConnectivityNode::'(.+)'"), item["ConductingEquipment.Terminals"][1]["Terminal.ConnectivityNode"]).captures[1]
                _recursive_dict_set!(data, [path..., i, "ConductingEquipment.BaseVoltage"], "BaseVoltage::'RAVENS_BaseV_$(base_voltages[cn]/1000.0)_kV'")
            end
        end
    end
end
