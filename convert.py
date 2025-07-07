import os
import argparse
import re
import math
import cantera as ct

R_CAL = 1.987

def format_value(val):
    return f"{val:.9g}"


def extract_true_third_body_contents(s: str):
    pattern = r"\+\s*(M)\b"
    
    returnValue = re.findall(pattern, s)
    # not empty, return the first match
    if returnValue:
        return returnValue[0]
    else:        
        return None
    
    
def extract_third_body_content(s: str):
    """
    Return a list of all the non-whitespace, non-')' substrings found inside
    parentheses in the input string.

    For example:
        extract_third_body_content("H + O2 (+AR) = HO2 (+AR)")
        # returns ['AR', 'AR']
    """
    pattern = r"\(\s*\+([^\)\s]+)\s*\)"
    
    returnValue = re.findall(pattern, s)
    # not empty, return the first match
    if returnValue:
        return returnValue[0]
    else:
        return None


def parse_transport(path):
    default = (1.0e-6, 120.0)
    table = {}
    if not path or not os.path.exists(path):
        return default, table
    current = None
    As = Ts = None
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('"') and line.endswith('"'):
                current = line.strip('"')
                As = Ts = None
            elif line.startswith('As'):
                As = float(line.split()[1].rstrip(';'))
            elif line.startswith('Ts'):
                Ts = float(line.split()[1].rstrip(';'))
            elif line.startswith('}') and current:
                if current == '.*':
                    default = (As, Ts)
                else:
                    table[current] = (As, Ts)
                current = None
    return default, table

def species_block(spec, default, table):
    thermo = spec.thermo
    coeffs = list(thermo.coeffs)
    Tcommon = coeffs[0]
    high = coeffs[1:8]
    low = coeffs[8:]
    # Cantera represents species which only provide a single NASA polynomial
    # using ``Tcommon`` equal to ``Thigh`` with duplicated coefficient blocks.
    # OpenFOAM, however, expects a reasonable midpoint temperature.  When
    # ``Tcommon`` coincides with ``Thigh`` and both polynomial blocks are
    # identical, fall back to a default value of ``1000 K`` which is commonly
    # used for NASA-7 fits.
    if (
        math.isclose(Tcommon, thermo.max_temp)
        and all(math.isclose(h, l, rel_tol=0, abs_tol=1e-12) for h, l in zip(high, low))
    ):
        Tcommon = 1000.0
    As, Ts = table.get(spec.name, default)
    lines = []
    lines.append(f"{spec.name}")
    lines.append("{")
    lines.append("    specie")
    lines.append("    {")
    lines.append(f"        molWeight       {format_value(spec.molecular_weight)};")
    lines.append("    }")
    lines.append("    thermodynamics")
    lines.append("    {")
    lines.append(f"        Tlow            {format_value(thermo.min_temp)};")
    lines.append(f"        Thigh           {format_value(thermo.max_temp)};")
    lines.append(f"        Tcommon         {format_value(Tcommon)};")
    lines.append("        highCpCoeffs    ( " + " ".join(format_value(v) for v in high) + " );")
    lines.append("        lowCpCoeffs     ( " + " ".join(format_value(v) for v in low) + " );")
    lines.append("    }")
    lines.append("    transport")
    lines.append("    {")
    lines.append(f"        As              {format_value(As)};")
    lines.append(f"        Ts              {format_value(Ts)};")
    lines.append("    }")
    lines.append("    elements")
    lines.append("    {")
    for el, amt in spec.composition.items():
        lines.append(f"        {el}               {format_value(amt)};")
    lines.append("    }")
    lines.append("}")
    return "\n".join(lines)


def writeSpecies(gas):
    """Return header lines listing all species in the mechanism."""
    species_names = gas.species_names
    lines = [
        "species         %d ( %s );" % (gas.n_species, " ".join(species_names)),
        "",
    ]
    return lines, species_names


def elements_block_lines(gas):
    """Return lines describing all elements used in the mechanism."""
    element_names = gas.element_names
    lines = ["elements", str(len(element_names)), "("]
    lines += element_names
    lines += [")", ";"]
    return lines


def writeThermo(gas, default_tp, table_tp, output_dir, species_names, header_lines=None):
    """Write the thermos file for all species.

    ``header_lines`` can be provided to prepend custom headers such as the
    species list.  If ``None``, no header is written.
    """
    thermo_lines = list(header_lines or [])
    for sp_name in species_names:
        sp = gas.species(sp_name)
        thermo_lines.append(species_block(sp, default_tp, table_tp))
        thermo_lines.append("")
    with open(os.path.join(output_dir, 'thermos'), 'w') as f:
        f.write("\n".join(thermo_lines) + "\n")


def writeReactions(
    gas,
    species_names,
    output_dir,
    pressure_atm=None,
    header_lines=None,
    elements_lines=None,
    rtype_suffix="",
):
    """Write the reactions file using all reactions from the mechanism.

    ``header_lines`` and ``elements_lines`` allow prepending species and
    elements information when required (eg, for OpenFOAM 7).
    """
    rxn_lines = []
    if header_lines:
        rxn_lines.extend(header_lines)
    if elements_lines:
        rxn_lines.extend(elements_lines)
        rxn_lines.append("")
    rxn_lines += ['reactions', '{']
    rxns = gas.reactions()
    i = 0
    out_idx = 0
    while i < len(rxns):
        rxn = rxns[i]
        try:
            if i + 1 < len(rxns) and is_reverse_pair(rxn, rxns[i + 1]):
                rxn_lines.append(
                    combined_reaction_block(
                        rxn,
                        rxns[i + 1],
                        out_idx,
                        species_names,
                        rtype_suffix=rtype_suffix,
                    )
                )
                i += 2
            else:
                rxn_lines.append(
                    reaction_block(
                        rxn,
                        out_idx,
                        species_names,
                        pressure_pa=pressure_atm * ct.one_atm if pressure_atm else None,
                        rtype_suffix=rtype_suffix,
                    )
                )
                i += 1
            out_idx += 1
        except NotImplementedError as e:
            print(f"Skipping reaction {i} ({rxn.equation}): {e}")
            i += 1
    rxn_lines.append('}')
    rxn_lines.append('Tlow            250;')
    rxn_lines.append('Thigh           5000;')
    with open(os.path.join(output_dir, 'reactions'), 'w') as f:
        f.write("\n".join(rxn_lines) + "\n")

def arrhenius_params(rate):
    A = rate.pre_exponential_factor
    beta = rate.temperature_exponent
    Ta = rate.activation_energy / ct.gas_constant
    return A, beta, Ta


def arrhenius_at_pressure(rates, pressure_pa):
    """Interpolate Arrhenius parameters for a given pressure.

    Parameters are interpolated in log space as implemented in ``test.py``.
    ``rates`` should be a list of ``(pressure, ArrheniusRate)`` tuples.
    """
    rates = sorted(rates, key=lambda pr: pr[0])

    for p, rate in rates:
        if abs(p - pressure_pa) / pressure_pa < 1e-6:
            return arrhenius_params(rate)

    for (p_i, rate_i), (p_j, rate_j) in zip(rates[:-1], rates[1:]):
        if p_i <= pressure_pa <= p_j:
            A_i, b_i, Ta_i = arrhenius_params(rate_i)
            A_j, b_j, Ta_j = arrhenius_params(rate_j)
            k = (math.log(pressure_pa) - math.log(p_i)) / (
                math.log(p_j) - math.log(p_i)
            )
            A = A_i * (A_j / A_i) ** k
            b = b_i + k * (b_j - b_i)
            Ta = Ta_i + k * (Ta_j - Ta_i)
            return A, b, Ta

    raise ValueError(
        "Pressure %.3g Pa outside range of provided PLOG data" % pressure_pa
    )

def format_equation(rxn):
    eq = rxn.equation.replace('<=>', '=').replace('=>', '=')
    eq = re.sub(r'\s*\(\+\s*[A-Za-z0-9_]+\s*\)', '', eq)
    eq = re.sub(r'\s*\+ M', '', eq)
    eq = re.sub(r'(\d+)\s+(\w)', r'\1\2', eq)
    eq = re.sub(r'\s+', ' ', eq).strip()

    orders = getattr(rxn, 'orders', {})
    if orders:
        pattern = re.compile(r'(\b\d*\.?\d*)([A-Za-z0-9_]+)')
        left, right = [s.strip() for s in eq.split('=')]

        def repl(match):
            coeff_str = match.group(1)
            species = match.group(2)
            stoich = float(coeff_str) if coeff_str else 1.0
            order = orders.get(species)
            if order is not None and abs(order - stoich) > 1e-9:
                return f"{coeff_str}{species}^{format_value(order)}"
            return f"{coeff_str}{species}"

        left = pattern.sub(repl, left)
        right = pattern.sub(repl, right)
        eq = f"{left} = {right}"

    return eq

def third_body_block(efficiencies, species_names, pure_third_body = True, must_return=False, only_one_third_body = False):
    # pure_third_body indicates third-body reactions without falloff and troe
    # one_third_body indicates only one of the species is considered in the third-body reactions

    if not efficiencies and not must_return:
        return []
    if not efficiencies:
        efficiencies = {sp: 1.0 for sp in species_names}

    default_value = 0.0 if only_one_third_body else 1.0
    

    lines = []
    if (not pure_third_body):
        lines += ["        thirdBodyEfficiencies",
                  "        {"]
    lines += [
        "            coeffs",
        "                " + str(len(species_names)),
        "            (",
    ]
    for sp in species_names:
        val = efficiencies.get(sp, default_value)
        lines.append(f"            ({sp} {format_value(val)})")
    lines += [
        "            )",
        "            ;",
    ]
    if (not pure_third_body):
        lines.append("        }")

    return lines

def k_block(name, rate):
    A, beta, Ta = arrhenius_params(rate)
    return [
        f"        {name}",
        "        {",
        f"            A               {format_value(A)};",
        f"            beta            {format_value(beta)};",
        f"            Ta              {format_value(Ta)};",
        "        }",
    ]

def k0_block(rate):
    return k_block("k0", rate)

def kinf_block(rate):
    return k_block("kInf", rate)

def F_block(alpha=None, T3=None, T1=None, T2=None):
    lines = ["        F", "        {"]
    if alpha is not None:
        lines.append(f"            alpha           {format_value(alpha)};")
    if T3 is not None:
        lines.append(f"            Tsss            {format_value(T3)};")
    if T1 is not None:
        lines.append(f"            Ts              {format_value(T1)};")
    if T2 is not None:
        lines.append(f"            Tss             {format_value(T2)};")
    lines.append("        }")
    return lines

def plog_block(rates):
    lines = ["        ArrheniusData", "        ("]
    for p, rate in rates:
        A, beta, Ta = arrhenius_params(rate)
        lines.append(
            f"            ({format_value(p)}  {format_value(A)} {format_value(beta)} {format_value(Ta)})"
        )
    lines += ["        )", "        ;"]
    return lines

def base_block(rtype, rxn, rate):
    A, beta, Ta = arrhenius_params(rate)
    return [
        f"        type            {rtype};",
        f"        reaction        \"{format_equation(rxn)}\";",
        f"        A               {format_value(A)};",
        f"        beta            {format_value(beta)};",
        f"        Ta              {format_value(Ta)};",
    ]

def reaction_block(
    rxn,
    index,
    species_names,
    pressure_pa=None,
    rtype_suffix="",
):
    only_one_third_body = is_one_third_body(rxn)
    
    prefix = 'reversible' if rxn.reversible else 'irreversible'
    if rxn.reaction_type == 'Arrhenius':
        rtype = f"{prefix}Arrhenius{rtype_suffix}"
        body = base_block(rtype, rxn, rxn.rate)
    elif rxn.reaction_type == 'three-body-Arrhenius':
        if rxn.input_data.get('type') == 'three-body':
            rtype = f"{prefix}ThirdBodyArrhenius{rtype_suffix}"
            body = base_block(rtype, rxn, rxn.rate)
            body += third_body_block(
                rxn.third_body.efficiencies,
                species_names,
                must_return=not rxn.third_body.efficiencies,
                only_one_third_body=only_one_third_body,
            )
        else:
            rtype = f"{prefix}Arrhenius{rtype_suffix}"
            body = base_block(rtype, rxn, rxn.rate)
    elif rxn.reaction_type == 'falloff-Troe':
        rtype = f"{prefix}ArrheniusTroeFallOff{rtype_suffix}"
        body = [f"        type            {rtype};", f"        reaction        \"{format_equation(rxn)}\";"]
        body += k0_block(rxn.rate.low_rate)
        body += kinf_block(rxn.rate.high_rate)

        if (len(rxn.rate.falloff_coeffs) == 3):
            alpha, T3, T1 = rxn.rate.falloff_coeffs
            body += F_block(alpha, T3, T1,4.5036e15)
        else:
            alpha, T3, T1, T2 = rxn.rate.falloff_coeffs
            body += F_block(alpha, T3, T1, T2)
        body += third_body_block(
            rxn.third_body.efficiencies,
            species_names,
            pure_third_body = False,
            must_return=not rxn.third_body.efficiencies,
            only_one_third_body=only_one_third_body,
        )
    elif rxn.reaction_type == 'falloff-Lindemann':
        rtype = f"{prefix}ArrheniusLindemannFallOff{rtype_suffix}"
        body = [f"        type            {rtype};", f"        reaction        \"{format_equation(rxn)}\";"]
        body += k0_block(rxn.rate.low_rate)
        body += kinf_block(rxn.rate.high_rate)
        body += F_block()
        body += third_body_block(
            rxn.third_body.efficiencies,
            species_names,
            pure_third_body = False,
            must_return=not rxn.third_body.efficiencies,
            only_one_third_body=only_one_third_body,
        )
    elif rxn.reaction_type == 'pressure-dependent-Arrhenius':
        if pressure_pa is None:
            rtype = f"{prefix}ArrheniusPLOG{rtype_suffix}"
            first_p, first_rate = rxn.rate.rates[0]
            body = base_block(rtype, rxn, first_rate)
            body += plog_block(rxn.rate.rates)
        else:
            rtype = f"{prefix}Arrhenius{rtype_suffix}"
            try:
                A, beta, Ta = arrhenius_at_pressure(rxn.rate.rates, pressure_pa)
                rate = ct.Arrhenius(A, beta, Ta * ct.gas_constant)
            except ValueError as exc:
                print(
                    f"Warning: {exc}. Using nearest tabulated pressure for {rxn.equation}."
                )
                # fallback to nearest pressure point
                rates = sorted(rxn.rate.rates, key=lambda pr: abs(pr[0] - pressure_pa))
                rate = rates[0][1]
            body = base_block(rtype, rxn, rate)
    else:
        raise NotImplementedError(f"Reaction type {rxn.reaction_type} not supported")

    lines = [f"    un-named-reaction-{index}", "    {"]
    lines.extend(body)
    lines.append("    }")
    return "\n".join(lines)

def is_reverse_pair(r1, r2):
    if r1.reaction_type != r2.reaction_type:
        return False
    if r1.reversible or r2.reversible:
        return False
    def same_stoich(a, b):
        if set(a.keys()) != set(b.keys()):
            return False
        for k in a:
            if abs(a[k] - b[k]) > 1e-12:
                return False
        return True
    return same_stoich(r1.reactants, r2.products) and same_stoich(r1.products, r2.reactants)

def combined_reaction_block(
    forward,
    reverse,
    index,
    species_names,
    rtype_suffix="",
):
    if forward.reaction_type == 'Arrhenius':
        rtype = f"nonEquilibriumReversibleArrhenius{rtype_suffix}"
    elif forward.reaction_type == 'three-body-Arrhenius':
        rtype = f"nonEquilibriumReversibleThirdBodyArrhenius{rtype_suffix}"
    else:
        raise NotImplementedError(f"Non-equilibrium reversible type for {forward.reaction_type} not supported")
    lines = [f"    un-named-reaction-{index}", "    {"]
    lines.append(f"        type            {rtype};")
    lines.append(f"        reaction        \"{format_equation(forward)}\";")
    lines.append("        forward")
    lines.append("        {")
    A, beta, Ta = arrhenius_params(forward.rate)
    lines.append(f"            A               {format_value(A)};")
    lines.append(f"            beta            {format_value(beta)};")
    lines.append(f"            Ta              {format_value(Ta)};")
    if forward.reaction_type == 'three-body-Arrhenius':
        lines += third_body_block(
            forward.third_body.efficiencies,
            species_names,
            must_return=not forward.third_body.efficiencies,
        )
    lines.append("        }")
    lines.append("        reverse")
    lines.append("        {")
    A, beta, Ta = arrhenius_params(reverse.rate)
    lines.append(f"            A               {format_value(A)};")
    lines.append(f"            beta            {format_value(beta)};")
    lines.append(f"            Ta              {format_value(Ta)};")
    if reverse.reaction_type == 'three-body-Arrhenius':
        lines += third_body_block(
            reverse.third_body.efficiencies,
            species_names,
            must_return=not reverse.third_body.efficiencies,
        )
    lines.append("        }")
    lines.append("    }")
    return "\n".join(lines)

def is_one_third_body(rxn):
    """
    Check if the reaction is a third-body reaction with only one species considered.
    """
    isLindemann = (rxn.reaction_type == 'falloff-Lindemann')
    isTroe = (rxn.reaction_type == 'falloff-Troe')
    isThirdBody = (rxn.reaction_type == 'three-body-Arrhenius')
    is_not_one_third_body = (extract_true_third_body_contents(rxn.equation) == "M")
    # if (is_not_one_third_body):
    #     print("\n")
    #     print("detect reaction:", rxn.equation)
    #     print(extract_true_third_body_contents(rxn.equation))
    if not (isLindemann or isTroe or isThirdBody):
        return False
    if is_not_one_third_body:
        return False
    if not rxn.reversible:
        return False
    else:
        third_body_species = extract_third_body_content(rxn.equation)
        if (third_body_species == "M"):
            return False
        else:
            # print("\n")
            # print(f"reaction type: {rxn.reaction_type}")
            # print(f"not M Detected third-body reaction: {rxn.equation}")
            # print(f"Third body species: {third_body_species}")
            return True


def default_converter(gas, default_tp, table_tp, output_dir, pressure_atm=None):
    """Converter logic for OpenFOAM 8 and newer."""
    header_lines, species_names = writeSpecies(gas)
    writeThermo(gas, default_tp, table_tp, output_dir, species_names, header_lines)
    with open(os.path.join(output_dir, 'elements'), 'w') as f:
        f.write("\n".join(elements_block_lines(gas)) + "\n")
    writeReactions(
        gas,
        species_names,
        output_dir,
        pressure_atm,
        rtype_suffix="",
    )


def of7_converter(gas, default_tp, table_tp, output_dir, pressure_atm=None):
    """Converter logic for OpenFOAM 7."""
    header_lines, species_names = writeSpecies(gas)
    # OF7 expects species and element information in the reactions file rather
    # than in the thermos file.
    writeThermo(gas, default_tp, table_tp, output_dir, species_names)
    elements_lines = elements_block_lines(gas)
    writeReactions(
        gas,
        species_names,
        output_dir,
        pressure_atm,
        header_lines=header_lines,
        elements_lines=elements_lines,
        rtype_suffix="Reaction",
    )


def convert(cantera_dir, output_dir, pressure_atm=None, version=10):
    os.makedirs(output_dir, exist_ok=True)
    mech = os.path.join(cantera_dir, 'chem.yaml')
    transport = os.path.join(cantera_dir, 'transportProperties')
    default_tp, table_tp = parse_transport(transport)
    gas = ct.Solution(mech)

    if int(version) == 7:
        of7_converter(gas, default_tp, table_tp, output_dir, pressure_atm)
    elif int(version) in (8, 9, 10):
        default_converter(gas, default_tp, table_tp, output_dir, pressure_atm)
    else:
        raise ValueError("Unsupported OpenFOAM version: %s" % version)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert Cantera YAML mechanism to OpenFOAM format')
    parser.add_argument('cantera_dir', help='directory containing chem.yaml and transportProperties')
    parser.add_argument('output_dir', help='directory where OpenFOAM files will be written')
    parser.add_argument(
        '--pressure',
        type=float,
        default=None,
        help='Pressure in atm for evaluating pressure-dependent reactions',
    )
    parser.add_argument(
        '--version',
        type=int,
        default=10,
        help='OpenFOAM major version (7-10)'
    )
    args = parser.parse_args()
    convert(
        args.cantera_dir,
        args.output_dir,
        pressure_atm=args.pressure,
        version=args.version,
    )
