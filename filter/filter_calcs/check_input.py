"""
Module to check the parsed input before starting the simulation.

Classes
=======

InputError:    does nothingm, parent: Exception

Functions
=========

check_input:    Checks the input parameters, logs errors and raises an
                InportError.
                Parameters:
                    parameters [struct]
@author: buechner_m <maria.buechner@gmail.com>
"""
import numpy as np
import parser_def
import materials
import logging
logger = logging.getLogger(__name__)

# %% Classes


class InputError(Exception):
    """
    InputError, parent 'Exception'
    """
    pass

# %% Public checking functions


def _test_check_parser(parameters):
    """
    Checking all input (from parser).

    Parameters
    ==========

    parameters [dict]

    """
    # Get parameter infos from parser, to link var_names and var_keys
    parser_info = parser_def.get_arguments_info(parser_def.input_parser())

    logger.info("Checking geometry input...")
    geometry_input(parameters, parser_info)
    logger.info("... done.")

    logger.info("Checking general input...")
    # % Minimal required input for all scenarios
    all_input(parameters, parser_info)
    logger.info("... done.")


def all_input(parameters, parser_info):
    """
    Checking all input necessary after geometries.

    Parameters
    ==========

    parameters [dict]
    parser_info [dict]:         parser_info[var_name] = [var_key, var_help]

    """
    try:
        # % Minimal required input for 'free', 'parallel', no gatings
        logger.debug("Checking general input...")

        # General input:
        logger.debug("Checking simulation input...")

        if not parameters['sampling_rate']:
            logger.debug("Sampling rate ({0}) is not specified, "
                         "set to pixel size * 1e-3."
                         .format(parser_info['sampling_rate'][0]))
            # Default to pixel_size *1e-3
            parameters['sampling_rate'] = parameters['pixel_size'] * 1e-3
            logger.debug("Sampling rate is {0} um, with pixel size {1} "
                         "um..".format(parameters['sampling_rate'],
                                       parameters['pixel_size']))
        if parameters['look_up_table'] == 'x0h' and \
                parameters['photo_only']:
            warning_message = ("With X0h material LUT cannot consider only "
                               "photo-absorption, resetting to 'False'.")
            logger.warn(warning_message)
            parameters['photo_only'] = False

        logger.debug("... done.")

        # Source:
        logger.debug("Checking Source input...")
        if parameters['beam_geometry'] == 'cone':
            if not parameters['focal_spot_size']:
                error_message = ("Input argument missing: 'focal_spot_size' "
                                 "({0})."
                                 .format(parser_info['focal_spot_size'][0]))
                logger.error(error_message)
                raise InputError(error_message)

        # Get spectrum
        [parameters['spectrum'], min_energy, max_energy] = \
            _get_spectrum(parameters['spectrum_file'],
                          parameters['spectrum_range'],
                          parameters['spectrum_step'],
                          parameters['design_energy'])  # [photons/pixel/sec]
        parameters['spectrum']['phontons'] = \
            parameters['spectrum']['phontons'] * \
            parameters['exposure_time'] #  [photons/pixel]
        # Material filter
        if parameters['thickness_filter'] and \
                not parameters['material_filter']:
            error_message = ("Filter material ({0}) must be specified."
                             .format(parser_info['material_filter'][0]))
            logger.error(error_message)
            raise InputError(error_message)
        if parameters['material_filter'] and \
                not parameters['thickness_filter']:
            error_message = ("Filter thickness ({0}) must be specified."
                             .format(parser_info['thickness_filter'][0]))
            logger.error(error_message)
            raise InputError(error_message)
        logger.debug("... done.")

        # Detector:
        logger.debug("Checking detector input...")
        # PSF right size?
        if parameters['detector_type'] == 'conv':
            if not parameters['point_spread_function']:
                error_message = ("Input argument missing: "
                                 "'point_spread_function' ({0})."
                                 .format(parser_info['point_spread_function']
                                         [0]))
                logger.error(error_message)
                raise InputError(error_message)
            if parameters['point_spread_function'] <= \
                    parameters['pixel_size']:
                # PSF too small, but be larger
                error_message = "PSF must be at larger than the pixel size."
                logger.error(error_message)
                raise InputError(error_message)

        # Threshold (error if > max energy and warninglog if < min)
        if parameters['detector_threshold']:
            if parameters['detector_threshold'] > max_energy:
                error_message = ("Detector threshold ({0}) must be <= the "
                                 "maximal energy ({1} keV)."
                                 .format(parser_info['detector_threshold'][0],
                                         max_energy))
                logger.error(error_message)
                raise InputError(error_message)
            elif parameters['detector_threshold'] < min_energy:
                logger.warning("Detector threshold ({0}) is smaller than the "
                               "minimal energy ){1} keV)."
                               .format(parser_info['detector_threshold'][0],
                                       min_energy))

        # Material thickness
        if parameters['thickness_detector'] and \
                not parameters['material_detector']:
            error_message = ("Detector material ({0}) must be specified."
                             .format(parser_info['material_detector'][0]))
            logger.error(error_message)
            raise InputError(error_message)
        if parameters['material_detector'] and \
                not parameters['thickness_detector']:
            error_message = ("Detector thickness ({0}) must be specified."
                             .format(parser_info['thickness_detector'][0]))
            logger.error(error_message)
            raise InputError(error_message)
        logger.debug("... done.")

        # Check all selected gratings (materials and phase/absorption)
        if 'G0' in parameters['component_list']:
            logger.debug("Checking G0...")
            _check_grating_input('g0', parameters, parser_info, False)
            logger.debug("... done.")
        if 'G1' in parameters['component_list']:
            logger.debug("Checking G1...")
            _check_grating_input('g1', parameters, parser_info, False)
            logger.debug("... done.")
        if 'G2' in parameters['component_list']:
            logger.debug("Checking G2...")
            _check_grating_input('g2', parameters, parser_info, False)
            logger.debug("... done.")

        # Check all materials if exist
        try:
            for var_name, value in parameters.iteritems():
                if 'material' in var_name and value:
                    materials.test_material(value, parameters['design_energy'],
                                            parameters['look_up_table'])
        except materials.MaterialError as e:
            error_message = ("Invalid material name in '{0}' ({1}): {2}"
                             .format(var_name,
                                     parser_info[var_name][0],
                                     str(e)))
            logger.error(error_message)
            raise InputError(error_message)

        logger.debug("... done.")  # General checking

    except AttributeError as e:
        error_message = "Input arguments missing: {}." \
                        .format(str(e).split()[-1])
        logger.error(error_message)
        raise InputError(error_message)


def geometry_input(parameters, parser_info):
    """
    Checking geometry input (everything to calculate the geometries).

    Parameters
    ==========

    parameters [dict]
    parser_info [dict]:         parser_info[var_name] = [var_key, var_help]

    """
    try:
        # % Minimal required input for 'free', 'parallel', no gatings
        logger.debug("Checking geometry input...")

        # GI Design:
        logger.debug("Checking GI input...")
        if parameters['gi_geometry'] != 'free':
            # If GI, talbot order necessary (unless dual_phase)
            if not parameters['talbot_order'] and not parameters['dual_phase']:
                error_message = ("Input argument missing: 'talbot_order' "
                                 "({0})."
                                 .format(parser_info['talbot_order'][0]))
                logger.error(error_message)
                raise InputError(error_message)
        else:
            if parameters['fixed_grating']:
                # If free geom and fixed_grating defined
                warning_message = ("Disregarding fixed grating choice for "
                                   "free geometry. Resetting to None.")
                logger.warning(warning_message)
                parameters['fixed_grating'] = None
        # Dual phase
        if parameters['dual_phase'] and parameters['gi_geometry'] != 'conv':
            error_message = ("Dual phase setup can only be calculated for "
                             "conventional geometry. Geometry is '{0}'."
                             .format(parameters['gi_geometry']))
            logger.error(error_message)
            raise InputError(error_message)
        # Store design wavelength
        parameters['design_wavelength'] = \
            materials.energy_to_wavelength(parameters['design_energy'])

        logger.debug("... done.")

        # Special scenarios:
        logger.debug("Checking geometry scenarios...")

        if (parameters['beam_geometry'] == 'parallel' and
                parameters['gi_geometry'] == 'sym') or \
                (parameters['beam_geometry'] == 'parallel' and
                 parameters['gi_geometry'] == 'inv'):
            error_message = ("Only '{0}' geometry valid for '{1}' beam "
                             "geometry."
                             .format(parameters['gi_geometry'],
                                     parameters['beam_geometry']))
            logger.error(error_message)
            raise InputError(error_message)

        parameters['component_list'] = ['Source', 'Detector']

        if parameters['beam_geometry'] == 'parallel':

            # Individual checks
            if parameters['gi_geometry'] == 'conv':
                # =============================================================
                # Conventional and parallel beam
                #
                # Conditions:
                #     no G0
                #     sample before or after G1 (bg1 or ag1)
                # Requirements:
                #     G1 and G2, one of them fixed
                # =============================================================
                logger.debug("Checking 'conv' geometry...")
                # Add G1 and G2
                parameters['component_list'].append('G1')
                parameters['component_list'].append('G2')
                # Sort updated component list
                parameters['component_list'].sort()
                # After sort, switch Source and Detector
                parameters['component_list'][0], \
                    parameters['component_list'][-1] = \
                    parameters['component_list'][-1], \
                    parameters['component_list'][0]

                # Fixed grating
                if parameters['fixed_grating'] == 'g0':
                    error_message = "The fixed grating must be either G1 or "
                    "G2."
                    logger.error(error_message)
                    raise InputError(error_message)
                elif not parameters['fixed_grating']:
                    error_message = ("Choose G1 or G2 as fixed grating ({0})."
                                     .format(parser_info['fixed_grating'][0]))
                    logger.error(error_message)
                    raise InputError(error_message)

                # Sample position (if defined)
                if parameters['sample_position']:
                    g1_index = parameters['component_list'].index('G1')
                    if parameters['sample_position'] == 'bg1':
                        parameters['component_list'].insert(g1_index,
                                                            'Sample')
                    elif parameters['sample_position'] == 'ag1':
                        parameters['component_list'].insert(g1_index+1,
                                                            'Sample')
                    else:
                        error_message = "Sample must be before or after G1."
                        logger.error(error_message)
                        raise InputError(error_message)
                logger.debug("... done.")
            else:
                # =============================================================
                # Free and parallel beam
                #
                # Conditions:
                #     no G0
                # Requirements:
                #     all distances between components must be given
                # =============================================================
                logger.debug("Checking 'free' geometry...")
                # Add all other components
                if parameters['type_g1']:
                    parameters['component_list'].append('G1')
                if parameters['type_g2']:
                    parameters['component_list'].append('G2')
                # Sort updated component list
                parameters['component_list'].sort()
                # After sort, switch Source and Detector
                parameters['component_list'][0], \
                    parameters['component_list'][-1] = \
                    parameters['component_list'][-1], \
                    parameters['component_list'][0]
                # Add sample
                if parameters['sample_position']:
                    if parameters['sample_position'] == 'as':
                        parameters['component_list'].insert(1, 'Sample')
                    elif parameters['sample_position'] == 'bg1':
                        reference_index = \
                            parameters['component_list'].index('G1')
                        parameters['component_list'].insert(reference_index,
                                                            'Sample')
                    elif parameters['sample_position'] == 'ag1':
                        reference_index = \
                            parameters['component_list'].index('G1')
                        parameters['component_list'].insert(reference_index+1,
                                                            'Sample')
                    elif parameters['sample_position'] == 'bg2':
                        reference_index = \
                            parameters['component_list'].index('G2')
                        parameters['component_list'].insert(reference_index,
                                                            'Sample')
                    elif parameters['sample_position'] == 'ag2':
                        reference_index = \
                            parameters['component_list'].index('G2')
                        parameters['component_list'].insert(reference_index+1,
                                                            'Sample')
                    else:
                        # 'bd'
                        parameters['component_list'].insert(-1, 'Sample')
                logger.debug("... done.")

            logger.debug("... done.")
        else:
            # Cone beam
            logger.debug("Checking 'cone' beam geometry...")
            if parameters['gi_geometry'] != 'free':
                logger.debug("Checking GI geometries...")
                # Common checks for not 'free' geometry
                # Add G1 and G2
                parameters['component_list'].append('G1')
                parameters['component_list'].append('G2')
                # G0
                if not parameters['type_g0']:
                    # No G0
                    if parameters['fixed_grating'] == 'g0':
                        if parameters['dual_pahse']:
                            error_message = "G0 is not defined, choose G1 "
                            "as fixed grating."
                        else:
                            error_message = "G0 is not defined, choose G1 or "
                            "G2 as fixed grating."
                        logger.error(error_message)
                        raise InputError(error_message)
                    elif (parameters['fixed_grating'] == 'g2' and
                          parameters['dual_pahse']):
                        error_message = "G1 must be fixed in dual phase "
                        "setup, not G2."
                        logger.error(error_message)
                        raise InputError(error_message)
                    elif not parameters['fixed_grating']:
                        if parameters['dual_pahse']:
                            error_message = ("Choose G1 as fixed grating "
                                             "({0})."
                                             .format(parser_info
                                                     ['fixed_grating'][0]))
                        else:
                            error_message = ("Choose G1 or G2 as fixed "
                                             "grating ({0})."
                                             .format(parser_info
                                                     ['fixed_grating'][0]))
                        logger.error(error_message)
                        raise InputError(error_message)
                    # Fixed distance
                    if parameters['gi_geometry'] != 'sym':
                        if not parameters['distance_source_g1'] and \
                                not parameters['distance_source_g2']:
                            error_message = ("Either distance from Source to "
                                             "G1 ({0}) OR Source to G2 ({1}) "
                                             "must be defined [mm]."
                                             .format(parser_info
                                                     ['distance_source_g1'][0],
                                                     parser_info
                                                     ['distance_source_g2']
                                                     [0]))
                            logger.error(error_message)
                            raise InputError(error_message)
                        elif parameters['distance_source_g1'] and \
                                parameters['distance_source_g2']:
                            if parameters['fixed_distance']:
                                logger.warning("Both distance from Source to "
                                               "G1 ({0}) AND Source to G2 "
                                               "({1}) are defined, choosing "
                                               "last choice of set distance "
                                               "({2})."
                                               .format(parser_info
                                                       ['distance_g0_g1'][0],
                                                       parser_info
                                                       ['distance_g0_g2'][0],
                                                       parameters
                                                       ['fixed_distance']))
                                fixed_distance = parameters['fixed_distance']
                            else:
                                # Fixed distance not defined
                                error_message = ("Either distance from Source "
                                                 "to G1 ({0}) OR Source to G2 "
                                                 "({1}) must be defined [mm]."
                                                 .format(parser_info
                                                         ['distance_g0_g1'][0],
                                                         parser_info
                                                         ['distance_g0_g2']
                                                         [0]))
                                logger.error(error_message)
                                raise InputError(error_message)
                        elif parameters['distance_source_g1']:
                            fixed_distance = 'distance_source_g1'
                        elif parameters['distance_source_g2']:
                            fixed_distance = 'distance_source_g2'
                    else:
                        fixed_distance = None
                else:
                    # With G0
                    # Add to component list (unless dual_phase)
                    if not parameters['dual_phase']:
                        parameters['component_list'].append('G0')
                    else:
                        error_message = ("Dual phase setup cannot include G0.")
                        logger.error(error_message)
                        raise InputError(error_message)
                    if not parameters['fixed_grating']:
                        error_message = ("Choose G0, G1 or G2 as fixed "
                                         "grating ({0})."
                                         .format(parser_info['fixed_grating']
                                                 [0]))
                        logger.error(error_message)
                        raise InputError(error_message)
                    # Fixed distance
                    if parameters['gi_geometry'] != 'sym':
                        if not parameters['distance_g0_g1'] and \
                                not parameters['distance_g0_g2']:
                            error_message = ("Either distance from G0 to G1 "
                                             "({0}) OR G0 to G2 ({1}) must be "
                                             "defined [mm]."
                                             .format(parser_info
                                                     ['distance_g0_g1'][0],
                                                     parser_info
                                                     ['distance_g0_g2'][0]))
                            logger.error(error_message)
                            raise InputError(error_message)
                        elif parameters['distance_g0_g1'] and \
                                parameters['distance_g0_g2']:
                            if parameters['fixed_distance']:
                                logger.warning("Both distance from G0 to G1 "
                                               "({0}) AND G0 to G2 ({1}) are "
                                               "defined, choosing last choice "
                                               "of set distance ({2})."
                                               .format(parser_info
                                                       ['distance_g0_g1'][0],
                                                       parser_info
                                                       ['distance_g0_g2'][0],
                                                       parameters
                                                       ['fixed_distance']))
                                fixed_distance = parameters['fixed_distance']
                            else:
                                # Fixed distance not defined
                                error_message = ("Either distance from G0 to "
                                                 "G1 ({0}) OR G0 to G2 ({1}) "
                                                 "must be defined [mm]."
                                                 .format(parser_info
                                                         ['distance_g0_g1'][0],
                                                         parser_info
                                                         ['distance_g0_g2']
                                                         [0]))
                                logger.error(error_message)
                                raise InputError(error_message)
                        elif parameters['distance_g0_g1']:
                            fixed_distance = 'distance_g0_g1'
                        elif parameters['distance_g0_g2']:
                            fixed_distance = 'distance_g0_g2'
                    else:
                        fixed_distance = None  # Sym
                parameters['fixed_distance'] = fixed_distance

                # Sort updated component list
                parameters['component_list'].sort()
                # After sort, switch Source and Detector
                parameters['component_list'][0], \
                    parameters['component_list'][-1] = \
                    parameters['component_list'][-1], \
                    parameters['component_list'][0]

                # Individaul checks
                if parameters['gi_geometry'] == 'conv':
                    # =========================================================
                    # Conventional and cone beam
                    #
                    # Conditions:
                    #     G0 optional
                    #     sample before G1 (bg1)
                    #     distance G2 to detector optional
                    #     distance from source (or G0) to G1
                    #       OR
                    #     distance from source (or G0) to G1
                    # Requirements:
                    #     G1 and G2
                    #     1 fixed grating
                    #
                    # =========================================================
                    logger.debug("Checking 'conv' geometry...")
                    # Sample position (if defined)
                    if parameters['sample_position']:
                        g1_index = parameters['component_list'].index('G1')
                        if parameters['sample_position'] == 'bg1':
                            parameters['component_list'].insert(g1_index,
                                                                'Sample')
                        else:
                            error_message = "Sample must be before G1."
                            logger.error(error_message)
                            raise InputError(error_message)
                    # Dual phase manual distances set?
                    if not parameters['distance_g1_g2']:
                        error_message = "Distance from G1 to G2 must be "
                        "defined."
                        logger.error(error_message)
                        raise InputError(error_message)
                    elif parameters['distance_g1_g2'] == 0:
                        error_message = "Distance from G1 to G2 must be "
                        "larger than 0."
                        logger.error(error_message)
                        raise InputError(error_message)
                    if not parameters['distance_g2_detector']:
                        error_message = "Distance from G2 to the detector "
                        "must be defined."
                        logger.error(error_message)
                        raise InputError(error_message)
                    elif parameters['distance_g2_detector'] == 0:
                        error_message = "Distance from G2 to the detector "
                        "must be larger than 0."
                        logger.error(error_message)
                        raise InputError(error_message)

                    logger.debug("... done.")
                elif parameters['gi_geometry'] == 'sym':
                    # =========================================================
                    # Symmetrical and cone beam
                    #
                    # Conditions:
                    #     G0 optional
                    #     sample before or after G1 (bg1 or ag1)
                    #     distance G2 to detector optional
                    #     distance from source (or G0) to G1
                    #       OR
                    #     distance from source (or G0) to G1
                    # Requirements:
                    #     G1 and G2
                    #     1 fixed grating
                    #
                    # =========================================================
                    logger.debug("Checking 'sym' geometry...")
                    # Sample position (if defined)
                    if parameters['sample_position']:
                        g1_index = parameters['component_list'].index('G1')
                        if parameters['sample_position'] == 'bg1':
                            parameters['component_list'].insert(g1_index,
                                                                'Sample')
                        elif parameters['sample_position'] == 'ag1':
                            parameters['component_list'].insert(g1_index+1,
                                                                'Sample')
                        else:
                            error_message = ("Sample must be before or after "
                                             "G1.")
                            logger.error(error_message)
                            raise InputError(error_message)
                    logger.debug("... done.")
                elif parameters['gi_geometry'] == 'inv':
                    # =========================================================
                    # Inverse and cone beam
                    #
                    # Conditions:
                    #     G0 optional
                    #     sample after G1 (ag1)
                    #     distance G2 to detector optional
                    #     distance from source (or G0) to G1
                    #       OR
                    #     distance from source (or G0) to G1
                    # Requirements:
                    #     G1 and G2
                    #     1 fixed grating
                    #
                    # =========================================================
                    logger.debug("Checking 'inv' geometry...")
                    # Sample position (if defined)
                    if parameters['sample_position']:
                        g1_index = parameters['component_list'].index('G1')
                        if parameters['sample_position'] == 'ag1':
                            parameters['component_list'].insert(g1_index+1,
                                                                'Sample')
                        else:
                            error_message = ("Sample must be before or after "
                                             "G1.")
                            logger.error(error_message)
                            raise InputError(error_message)
                    logger.debug("... done.")
                logger.debug("... done.")
            else:
                # =============================================================
                # Free and cone beam
                #
                # Conditions:
                #
                # Requirements:
                #     all distances between components must be given
                # =============================================================
                logger.debug("Checking 'free' geometry...")
                # Add all other components
                if parameters['type_g0']:
                    parameters['component_list'].append('G0')
                if parameters['type_g1']:
                    parameters['component_list'].append('G1')
                if parameters['type_g2']:
                    parameters['component_list'].append('G2')
                # Sort updated component list
                parameters['component_list'].sort()
                # After sort, switch Source and Detector
                parameters['component_list'][0], \
                    parameters['component_list'][-1] = \
                    parameters['component_list'][-1], \
                    parameters['component_list'][0]
                # Add sample
                if parameters['sample_position']:
                    if parameters['sample_position'] == 'as':
                        parameters['component_list'].insert(1, 'Sample')
                    elif parameters['sample_position'] == 'bg0':
                        reference_index = \
                            parameters['component_list'].index('G0')
                        parameters['component_list'].insert(reference_index,
                                                            'Sample')
                    elif parameters['sample_position'] == 'ag0':
                        reference_index = \
                            parameters['component_list'].index('G0')
                        parameters['component_list'].insert(reference_index+1,
                                                            'Sample')
                    elif parameters['sample_position'] == 'bg1':
                        reference_index = \
                            parameters['component_list'].index('G1')
                        parameters['component_list'].insert(reference_index,
                                                            'Sample')
                    elif parameters['sample_position'] == 'ag1':
                        reference_index = \
                            parameters['component_list'].index('G1')
                        parameters['component_list'].insert(reference_index+1,
                                                            'Sample')
                    elif parameters['sample_position'] == 'bg2':
                        reference_index = \
                            parameters['component_list'].index('G2')
                        parameters['component_list'].insert(reference_index,
                                                            'Sample')
                    elif parameters['sample_position'] == 'ag2':
                        reference_index = \
                            parameters['component_list'].index('G2')
                        parameters['component_list'].insert(reference_index+1,
                                                            'Sample')
                    else:
                        # 'bd'
                        parameters['component_list'].insert(-1, 'Sample')
                logger.debug("... done.")

            logger.debug("... done.")

        # Set optional distances from None to 0
        if parameters['beam_geometry'] == 'cone' and \
                parameters['gi_geometry'] != 'free':
            # G2 to detector
            if parameters['distance_g2_detector'] is None:
                logger.debug("Setting undefined optional distance "
                             "'distance_g2_detector to: 0")
                parameters['distance_g2_detector'] = 0.0

            if 'G0' in parameters['component_list'] and \
                    parameters['distance_source_g0'] is None:
                logger.debug("Setting undefined optional distance "
                             "'distance_source_g0' to: 0")
                parameters['distance_source_g0'] = 0.0

        # Info
        logger.info("Beam geometry is '{0}' and setup geometry is '{1}'."
                    .format(parameters['beam_geometry'],
                            parameters['gi_geometry']))
        logger.info("Setup consists of: {0}."
                    .format(parameters['component_list']))
        if 'Sample' not in parameters['component_list']:
            logger.info("No sample included.")

        # Check fixed grating
        if parameters['gi_geometry'] != 'free':
            logger.info("Fixed grating is: '{0}'."
                        .format(parameters['fixed_grating']))
            # Fixed grating
            logger.debug("Checking {0}...".format(parameters['fixed_grating']))
            _check_grating_input(parameters['fixed_grating'], parameters,
                                 parser_info, True)
            logger.debug("... done.")
            # Check G2 if dual phase as if it was fixed for phase shift
            # settings (G1 is fixed)
            if parameters['dual_phase']:
                logger.debug("Checking G2...")
                fixed_grating = parameters['fixed_grating']
                parameters['fixed_grating'] = 'g2'
                _check_grating_input('g2', parameters, parser_info, True)
                parameters['fixed_grating'] = fixed_grating
                logger.debug("... done.")
            # Fixed distance
            if parameters['beam_geometry'] == 'cone':
                logger.info("Fixed distance is: {0}."
                            .format(fixed_distance))

        # Check remaining components
        # Sample distance, shape, material etc.
        logger.debug("Checking remaining components...")
        if 'Sample' in parameters['component_list']:
            logger.debug("Checking sample input...")
            if not parameters['sample_distance']:
                error_message = ("Distance from sample to reference component "
                                 "must be specified.")
                logger.error(error_message)
                raise InputError(error_message)
            # ########## Temp ########################
            if not parameters['sample_diameter']:
                error_message = "Sample diameter must be specified."
                logger.error(error_message)
                raise InputError(error_message)
            # ########## Temp ########################
            logger.debug("... done.")

        if parameters['gi_geometry'] == 'free':
            # Chack all necessary distances
            logger.debug("Checking distances for 'free' input...")
            for index, component in \
                    enumerate(parameters['component_list'][:-1]):
                current_distance = ('distance_' + component.lower() + '_' +
                                    parameters['component_list'][index+1]
                                    .lower())
                if not parameters[current_distance]:
                    error_message = ("{0} ({1}) not defined."
                                     .format(parser_info[current_distance][1]
                                             .split('.')[0],
                                             parser_info[current_distance][0]))
                    logger.error(error_message)
                    raise InputError(error_message)

            logger.debug("... done.")

        # Component checking done
        logger.debug("... done.")
        # Scenarios done
        logger.debug("... done.")
        # General checking
        logger.debug("... done.")

    except AttributeError as e:
        error_message = "Input arguments missing: {}." \
                        .format(str(e).split()[-1])
        logger.error(error_message)
        raise InputError(error_message)

# %% Public utility functions


# %% Private checking functions


# %% Private utility functions

def _get_spectrum(spectrum_file, range_, spectrum_step, design_energy):
    """
    Load spectrum from file or define based on range (min, max). Returns
    energies and relative photons (photons/pixel/sec).

    Parameters
    ==========

    spectrum_file:              path to spectrum file
    range_ [keV, keV]:          [min, max]
    spectrum_step [keV]
    design_energy [keV]

    Returns
    =======

    [spectrum, min, max]        spectrum:   [energies, photons]
                                            [keV, photons/pixel/sec]
                                min: minimal energy [keV]
                                max: maximal energy [keV]

    Notes
    =====

    if from file and range:
        range is set within loaded spectrum. Photons not rescaled.
    if range:
        range from min to max, step 1 keV. Homogenuous photons distribution.

    """
    # Read from file
    if spectrum_file is not None:
        logger.info("Reading spectrum from file at:\n{}..."
                    .format(spectrum_file))
        spectrum = _read_spectrum(spectrum_file)
        # Set range
        if range_ is not None:
            # Min and max in right order?
            if range_[1] <= range_[0]:
                error_message = ("Energy range maximum value ({0} keV) must "
                                 "be larger than minimum value ({1} keV)."
                                 .format(range_[0], range_[1]))
                logger.error(error_message)
                raise InputError(error_message)
            # Check if within bounds of spectrum
            if range_[0] >= max(spectrum['energies']):
                error_message = ("Energy range minimum value must be smaller "
                                 "than spectrum maximum ({0} keV)."
                                 .format(max(spectrum['energies'])))
                logger.error(error_message)
                raise InputError(error_message)
            if range_[1] <= min(spectrum['energies']):
                error_message = ("Energy range maximum value must be larger "
                                 "than spectrum minimum ({0} keV)."
                                 .format(min(spectrum['energies'])))
                logger.error(error_message)
                raise InputError(error_message)

            # Find min and max closes to range min and max
            [min_energy, min_index] = _nearest_value(spectrum['energies'],
                                                     range_[0])
            [max_energy, max_index] = _nearest_value(spectrum['energies'],
                                                     range_[1])
            # More than 1 energy in spectrum?
            if min_energy == max_energy:
                error_message = ("Energy minimum value same as maximum. Range"
                                 "minimum {0} and maximum {1} too close "
                                 "together.".format(range_[0], range_[1]))
                logger.error(error_message)
                raise InputError(error_message)
            spectrum['energies'] = spectrum['energies'][min_index:max_index+1]
            spectrum['photons'] = spectrum['photons'][min_index:max_index+1]
            logger.debug("\tSet energy range from {0} to {1} keV."
                         .format(min_energy, max_energy))
        logger.info("... done.")
    # Check range input
    elif range_ is not None:
        # Min and max in right order?
        if range_[1] <= range_[0]+spectrum_step:
            error_message = ("Energy range maximum value ({0} keV) must be "
                             "at least {1} keV larger than minimum value "
                             "({2} keV)."
                             .format(range_[0], spectrum_step, range_[1]))
            logger.error(error_message)
            raise InputError(error_message)
        # Calc spectrum
        spectrum = dict()
        # Calc from range
        logger.info("Setting spectrum based on input range...")
        spectrum['energies'] = np.arange(range_[0],
                                         range_[1]+spectrum_step,
                                         spectrum_step,
                                         dtype=np.float)
        spectrum['photons'] = np.ones(len(spectrum['energies']),
                                      dtype=np.float)
        logger.debug("\tSet all photons to {}."
                     .format(spectrum['photons'][0]))
        # Convert to struct
        logger.info("... done.")
    # Both spectrum_file and _range are None, use design energy as spectrum
    else:
        logger.info("Only design energy specified, calculating only for {} "
                    "keV...".format(design_energy))
        spectrum = dict()
        spectrum['energies'] = np.array(design_energy, dtype=np.float)
        spectrum['photons'] = np.array(1, dtype=np.float)
        logger.debug("\tSet photons to 1.")
        logger.info("... done.")
        logger.info("Spectrum is design energy {} keV."
                    .format(spectrum['energies']))
        return spectrum, spectrum['energies'], spectrum['energies']

    # Check and show spectrum results
    min_energy = min(spectrum['energies'])
    max_energy = max(spectrum['energies'])
    # Design energy in spectrum?
    if design_energy < min_energy or \
       design_energy > max_energy:
        error_message = ("Design energy ({0} keV) must be within "
                         "spectrum range (min: {1} kev, max: {2} keV).") \
                         .format(design_energy,
                                 min_energy, max_energy)
        logger.error(error_message)
        raise InputError(error_message)
    logger.debug("Design energy within spectrum.")
    logger.info("Spectrum from {0} keV to {1} keV in {2} keV steps."
                .format(min_energy, max_energy,
                        spectrum['energies'][1]-spectrum['energies'][0]))
    return spectrum, min_energy, max_energy


def _read_spectrum(spectrum_file_path):
    """
    Read from spectrum file.

    Parameters
    ==========

    spectrum_file_path [str]:       to .csv or .txt file.
                                    Delimiter is ',', see Format

    Returns
    =======

    spectrum [dict] [keV]:          spectrum['energies']
                                    spectrum['photons']

    Format
    ======

    energy, photons
    1, 10e3
    2, 23e4
    .,.
    .,.
    .,.

    """
    # Read dict from file
    logger.debug("Reading from file {}...".format(spectrum_file_path))
    spectrum_struct_array = np.genfromtxt(spectrum_file_path, delimiter=',',
                                          names=True)  # np ndarray
    if 'energy' in spectrum_struct_array.dtype.names:
        # Rename 'energy' to 'energies'
        spectrum_struct_array.dtype.names = ('energies', 'photons')
    # Convert to dict
    spectrum = dict()
    try:
        spectrum['energies'] = spectrum_struct_array['energies']
        spectrum['photons'] = spectrum_struct_array['photons']
    except AttributeError as e:
        error_message = "Spectrum file at {0} is missing '{1}'-column." \
                        .format(spectrum_file_path, str(e).split()[-1])
        logger.error(error_message)
        raise InputError(error_message)

    # Check if more than 2 energies in spectrum
    if len(spectrum['energies']) <= 1:
                error_message = ("Spectrum file only contains 1 energy."
                                 "Minimum is 2.")
                logger.error(error_message)
                raise InputError(error_message)
    logger.debug("... done.")
    return spectrum


def _nearest_value(array, value):
    """
    Funtion to find the nearest value of a number within a numpy array.

    Parameters
    ==========

    array [numpy array]     array to be searched
    value                   target number

    Returns
    =======

    [nearest_value, index]

    """
    nearest_index = (np.abs(array-value)).argmin()
    return array[nearest_index], nearest_index


def _check_grating_input(grating, parameters, parser_info, geometry):
    """
    Check grating input.

    Parameters
    ==========

    grating [str]:          converts to lower case, 'g0', 'g1' or 'g2'
    parameters [dict]
    parser_info [dict]:     parser_info[var_name] = [var_key, var_help]
    geometry [boolean]:     called from check geometry (True), or check all (F)

    Notes
    =====

    Fixed grating is always None for free geometry, see
    check_input.geometry_input.

    Basic required input:
        type (is defined)
        pitch (if fixed grating or free geoemtry, except G2 and dual phase)
        duty cycle (if fixed grating or free geometry, except G2 and dual
                    phase)

        material (optional for geometry calc)
        Phase and/or thickness and/or absorption: (cals only if not geometry)
            If free:
                abs:
                    required: thickness or absorption (dominant)
                    calculated: missing thickness or absorption; phase=None
                phase:
                    required: phase (dominant) or thickness
                    calculated: missing phase or thickness; absorption=None
                mix:
                    required: thickness (dominant) or phase or absorption
                    calculated: missing two thickness or phase or absorption

            If GI:
                if G0/G2: (if dual phase: G2: phase or mix, like G1)
                    abs:
                        required: thickness or absorption (dominant)
                        calculated: missing thickness or absoprtion; phase=None
                    mix:
                        required: thickness (dominant) or absorption
                        calculated: missing thickness or absorption and phase
                if G1:
                    phase:
                        required: phase
                        calculated: thickness; absorption=None
                    mix:
                        required: phase
                        calculated: thickness and absorption

            (-> abs same for free and GI)

    """
    grating = grating.lower()
    # Is defined?
    if not parameters['type_'+grating]:
        error_message = ("Type of {0} ({1}) not defined."
                         .format(grating.upper(),
                                 parser_info['type_'+grating][0]))
        logger.error(error_message)
        raise InputError(error_message)

    # Check grating types for GI setups (optional for geometry calc)
    if parameters['gi_geometry'] != 'free' and not geometry:
        # G0 (abs or mix)
        if grating == 'g0' and parameters['type_'+grating] == 'phase':
            error_message = ("Type of G0 ({0}) must be 'mix' or 'abs'."
                             .format(parser_info['type_'+grating][0]))
            logger.error(error_message)
            raise InputError(error_message)
        # G1 (phase or mix)
        if grating == 'g1' and parameters['type_'+grating] == 'abs':
            error_message = ("Type of G1 ({0}) must be 'mix' or 'phase'."
                             .format(parser_info['type_'+grating][0]))
            logger.error(error_message)
            raise InputError(error_message)
        # G2 (abs or mix for classic GI, phase or mix for dual phase)
        if grating == 'g2' and not parameters['dual_phase'] and \
                parameters['type_'+grating] == 'phase':
            error_message = ("Type of G2 ({0}) must be 'mix' or 'abs'."
                             .format(parser_info['type_'+grating][0]))
            logger.error(error_message)
            raise InputError(error_message)
        elif grating == 'g2' and parameters['dual_phase'] and \
                parameters['type_'+grating] == 'abs':
            error_message = ("Type of G2 ({0}) must be 'mix' or 'phase'."
                             .format(parser_info['type_'+grating][0]))
            logger.error(error_message)
            raise InputError(error_message)

    # Basic required input (except G2 and dual phase)
    # If fixed grating (for none-free input) or free input
    if ((parameters['gi_geometry'] != 'free' and
            grating == parameters['fixed_grating']) or
            parameters['gi_geometry'] == 'free') and not \
            (parameters['dual_phase'] and grating == 'g2'):
        if not parameters['pitch_'+grating]:
            error_message = ("Pitch of {0} ({1}) must be defined."
                             .format(grating.upper(),
                                     parser_info['pitch_'+grating][0]))
            logger.error(error_message)
            raise InputError(error_message)
        if not parameters['duty_cycle_'+grating]:
            error_message = ("Duty cycle of {0} ({1}) must be defined."
                             .format(grating.upper(),
                                     parser_info['duty_cycle_'+grating][0]))
            logger.error(error_message)
            raise InputError(error_message)
        elif parameters['duty_cycle_'+grating] <= 0 or \
                parameters['duty_cycle_'+grating] >= 1:
            error_message = ("Duty cycle of {0} ({1}) must be within ]0...1[."
                             .format(grating.upper(),
                                     parser_info['duty_cycle_'+grating][0]))
            logger.error(error_message)
            raise InputError(error_message)

    # Phase and/or thickness and/or absorption
    if parameters['type_'+grating] == 'abs':
        # Absorption grating
        if not parameters['thickness_'+grating]:
            error_message = ("Thickness of {0} ({1}) must be defined."
                             .format(grating.upper(),
                                     parser_info['thickness_'+grating][0]))
            logger.error(error_message)
            raise InputError(error_message)
        if parameters['phase_shift_'+grating]:
            warning_message = ("Phase shift of {0} is defined, but ignored."
                               .format(grating.upper()))
            logger.warn(warning_message)
            parameters['phase_shift_'+grating] = None
    elif parameters['type_'+grating] == 'phase':
        # Phase grating
        if parameters['gi_geometry'] == 'free':
            # Free input
            if not parameters['thickness_'+grating] and \
                    not parameters['phase_shift_'+grating]:
                # Nothing defined
                error_message = ("Thickness ({0}) OR phase shift ({1}) of {2} "
                                 "must be defined."
                                 .format(parser_info['thickness_'+grating][0],
                                         parser_info['phase_shift_'+grating]
                                         [0],
                                         grating.upper()))
                logger.error(error_message)
                raise InputError(error_message)
            if parameters['phase_shift_'+grating]:
                # Phase defined
                if parameters['thickness_'+grating]:
                    # Thickness as well
                    warning_message = ("Thickness AND phase shift of {0} are "
                                       "defined. Basing calculations on phase "
                                       "shift.".format(grating.upper()))
                    logger.warn(warning_message)
                # Calc thickness
                if not geometry:
                    parameters['thickness_'+grating] = \
                        materials.shift_to_height(parameters
                                                  ['phase_shift_'+grating],
                                                  parameters
                                                  ['material_'+grating],
                                                  parameters['design_energy'],
                                                  photo_only=
                                                  parameters['photo_only'],
                                                  source=
                                                  parameters['look_up_table'])
            else:
                # Phase not defined, but thickness
                # Calc phase shift
                if not geometry:
                    parameters['phase_shift_'+grating] = \
                        materials.height_to_shift(parameters
                                                  ['thickness_'+grating],
                                                  parameters
                                                  ['material_'+grating],
                                                  parameters['design_energy'],
                                                  photo_only=
                                                  parameters['photo_only'],
                                                  source=
                                                  parameters['look_up_table'])
        else:
            # GI setup
            # Either phase G1 (normal and dual phase) or
            # phase G2 (for dual phase)
            if not parameters['phase_shift_'+grating]:
                if parameters['dual_phase']:
                    error_message = ("Phase shift ({0}) of {1} "
                                     "must be defined."
                                     .format(parser_info['phase_shift_' +
                                                         grating][0],
                                             grating.upper()))
                else:
                    error_message = ("Phase shift ({0}) of {1} "
                                     "must be defined as pi or pi/2."
                                     .format(parser_info['phase_shift_' +
                                                         grating][0],
                                             grating.upper()))
                logger.error(error_message)
                raise InputError(error_message)

            if (not parameters['dual_phase'] and
                not((round(parameters['phase_shift_'+grating] - np.pi) == 0) or
                   (round(parameters['phase_shift_'+grating] - np.pi/2) == 0))):
                error_message = ("Phase shift ({0}) of {1} must be 'pi' or "
                                 "'pi/2'."
                                 .format(parser_info['phase_shift_'+grating]
                                         [0],
                                         grating.upper()))
                logger.error(error_message)
                raise InputError(error_message)

            if parameters['thickness_'+grating] and \
                    parameters['phase_shift_'+grating]:
                warning_message = ("Thickness AND phase shift of {0} are "
                                   "defined. Basing calculations on phase "
                                   "shift.".format(grating.upper()))
                logger.warn(warning_message)
            # Calc thickness
            if not geometry:
                parameters['thickness_'+grating] = \
                    materials.shift_to_height(parameters['phase_shift_' +
                                                         grating],
                                              parameters['material_'+grating],
                                              parameters['design_energy'],
                                              photo_only=
                                              parameters['photo_only'],
                                              source=
                                              parameters['look_up_table'])
    else:
        # Mix grating
        if parameters['gi_geometry'] == 'free':
            # Free input
            if not parameters['thickness_'+grating] and \
                    not parameters['phase_shift_'+grating]:
                # Nothing defined
                error_message = ("Thickness ({0}) OR phase shift ({1}) of {2} "
                                 "must be defined."
                                 .format(parser_info['thickness_'+grating][0],
                                         parser_info['phase_shift_'+grating]
                                         [0],
                                         grating.upper()))
                logger.error(error_message)
                raise InputError(error_message)
            if parameters['thickness_'+grating]:
                # Thickness defined
                if parameters['phase_shift_'+grating]:
                    # Phase as well
                    warning_message = ("Thickness AND phase shift of {0} are "
                                       "defined. Basing calculations on "
                                       "thickness.".format(grating.upper()))
                    logger.warn(warning_message)
                # Calc phase shift
                if not geometry:
                    parameters['phase_shift_'+grating] = \
                        materials.height_to_shift(parameters
                                                  ['thickness_'+grating],
                                                  parameters
                                                  ['material_'+grating],
                                                  parameters['design_energy'],
                                                  photo_only=
                                                  parameters['photo_only'],
                                                  source=
                                                  parameters['look_up_table'])
            else:
                # Phase shift defined
                # Calc thickness
                if not geometry:
                    parameters['thickness_'+grating] = \
                        materials.shift_to_height(parameters
                                                  ['phase_shift_'+grating],
                                                  parameters['material_' +
                                                             grating],
                                                  parameters['design_energy'],
                                                  photo_only=
                                                  parameters['photo_only'],
                                                  source=
                                                  parameters['look_up_table'])
        else:
            # GI setup
            if grating == 'g0' or (grating == 'g2' and
                                   not parameters['dual_phase']):
                # For G0 or G2 if not dual_phase
                if not parameters['thickness_'+grating]:
                    # Thickness not defined
                    error_message = ("Thickness ({0}) of {1} "
                                     "must be defined."
                                     .format(parser_info['thickness_'+grating]
                                             [0],
                                             grating.upper()))
                    logger.error(error_message)
                    raise InputError(error_message)
                if parameters['phase_shift_'+grating]:
                    # Phase as well
                    warning_message = ("Thickness AND phase shift of {0} are "
                                       "defined. Basing calculations on "
                                       "thickness.".format(grating.upper()))
                    logger.warn(warning_message)
                # Calc phase shift
                if not geometry:
                    parameters['phase_shift_'+grating] = \
                        materials.height_to_shift(parameters
                                                  ['thickness_'+grating],
                                                  parameters
                                                  ['material_'+grating],
                                                  parameters['design_energy'],
                                                  photo_only=
                                                  parameters['photo_only'],
                                                  source=
                                                  parameters['look_up_table'])
            else:
                # G1 (normal and dual phase) or G2 if dual phase
                if not parameters['phase_shift_'+grating]:
                    error_message = ("Phase shift ({0}) of {1} "
                                     "must be defined."
                                     .format(parser_info
                                             ['phase_shift_'+grating][0],
                                             grating.upper()))
                    logger.error(error_message)
                    raise InputError(error_message)
                if parameters['phase_shift_'+grating]:
                    # Phase as well
                    warning_message = ("Thickness AND phase shift of {0} are "
                                       "defined. Basing calculations on "
                                       "thickness.".format(grating.upper()))
                    logger.warn(warning_message)
                # Calc thickness
                if not geometry:
                    parameters['thickness_'+grating] = \
                        materials.shift_to_height(parameters
                                                  ['phase_shift_'+grating],
                                                  parameters['material_' +
                                                             grating],
                                                  parameters['design_energy'],
                                                  photo_only=
                                                  parameters['photo_only'],
                                                  source=
                                                  parameters['look_up_table'])

    # Optional input if not geometry check
    # Always required
    if not parameters['material_'+grating] and not geometry:
        error_message = ("Material of {0} ({1}) must be defined."
                         .format(grating.upper(),
                                 parser_info['material_'+grating][0]))
        logger.error(error_message)
        raise InputError(error_message)

    # Optional input
    # Wafer
    if parameters['wafer_thickness_'+grating] and \
            not parameters['wafer_material_'+grating]:
        error_message = ("Wafer material of {0} ({1}) must be specified."
                         .format(grating.upper(),
                                 parser_info['wafer_material_'+grating][0]))
        logger.error(error_message)
        raise InputError(error_message)
    if parameters['wafer_material_'+grating] and \
            not parameters['wafer_thickness_'+grating]:
        error_message = ("Wafer thickness of {0} ({1}) must be specified."
                         .format(grating.upper(),
                                 parser_info['wafer_thickness_'+grating][0]))
        logger.error(error_message)
        raise InputError(error_message)
    # Grating fill
    if parameters['fill_thickness_'+grating] and \
            not parameters['fill_material_'+grating]:
        error_message = ("Fill material of {0} ({1}) must be specified."
                         .format(grating.upper(),
                                 parser_info['fill_material_'+grating][0]))
        logger.error(error_message)
        raise InputError(error_message)
    if parameters['fill_material_'+grating] and \
            not parameters['fill_thickness_'+grating]:
        error_message = ("Fill thickness of {0} ({1}) must be specified."
                         .format(grating.upper(),
                                 parser_info['fill_thickness_'+grating][0]))
        logger.error(error_message)
        raise InputError(error_message)

    # Shape of grating
    if parameters[grating+'_bent']:
        if parameters['beam_geometry'] == 'parallel':
            warning_message = ("{0} is bent in a parallel beam geometry. "
                               "Ignoring bent grating parameters..."
                               .format(grating.upper()))
            logger.warning(warning_message)
            parameters[grating+'_bent'] = False
            parameters['radius_'+grating] = None
            parameters[grating+'_matching'] = False
        elif parameters[grating+'_matching'] and parameters['radius_'+grating]:
            warning_message = ("{0} is bent with matching radius AND a radius "
                               "is set. Ignoring set radius..."
                               .format(grating.upper()))
            logger.warning(warning_message)
            parameters['radius_'+grating] = None
        elif not parameters[grating+'_matching'] and \
                not parameters['radius_'+grating]:
            error_message = ("Radius of bent {0} is required."
                             .format(grating.upper()))
            logger.error(error_message)
            raise InputError(error_message)
    else:
        if parameters[grating+'_matching']:
            warning_message = ("{0} is straight and cannot match its distance "
                               "from source. Ignoring matching flag..."
                               .format(grating.upper()))
            logger.warning(warning_message)
            parameters[grating+'_matching'] = False
        if parameters['radius_'+grating]:
            warning_message = ("{0} is straight and does not have a radius. "
                               "gnoring set radius..."
                               .format(grating.upper()))
            logger.warning(warning_message)
            parameters['radius_'+grating] = None