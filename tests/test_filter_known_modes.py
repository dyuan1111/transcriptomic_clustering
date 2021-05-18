import pytest

import numpy as np
import pandas as pd
import scanpy as sc

import transcriptomic_clustering as tc

@pytest.fixture
def genes():
    """
        genes
    """
    genes = ['Plp1', 'Npy', 'Cnp', 'Mal', 'Trf', 'Enpp2', 'Penk', 'Cnr1', 'Cd9',
                'Rgs5', 'Krt73', 'Myl1', 'Lpar1', 'Arhgap36', 'Vip', 'Rgs4',
                'Mobp', 'Tac2', '6330527O06Rik', 'Dmp1']
    
    return genes

@pytest.fixture
def principal_components(genes):
    """
        principal components
    """
    mat_pcs = np.array([[0.002857591636481942, -0.0354102055296144, 0.03249323043004619,
        0.055205398305197745, 0.04305623050772404],
       [-0.4240749068725877, 0.3955913078965465, -0.05559982402976456,
        0.34183306838972055, -0.4306519684095824],
       [-0.010365029305356483, -0.10828103360262205,
        0.021531245294265332, 0.04029950841330414, -0.04690326670059291],
       [0.0008139138770373759, 0.0015875831615960966,
        0.02828577591205593, -0.027742334187004246,
        -0.023663637678693082],
       [-0.0003600030414637642, -0.03240472172146662,
        0.004780957009875549, 0.003083358753351117,
        -0.004709799560498136],
       [-0.04171670245104074, -0.054264449676871426, 0.3989904532820567,
        0.25359079868413875, -0.7159660256802589],
       [0.39982611321528977, -0.12764441080640382, -0.14143664558942642,
        -0.023963648540285133, 0.4573350489202506],
       [0.05378849672013527, 0.03914426964853044, 0.3092658322561669,
        -0.5048850434413988, 0.12142840351519307],
       [-0.007169182817480407, -0.03405309463464166, 0.04666958798795709,
        0.01736164494235163, 0.0055572088279958395],
       [-0.019864084169154535, 0.0058319878150946145,
        -0.04226176748153477, -0.04468406823299311,
        -0.039769620541415786],
       [-0.025394541899113793, -0.06489107975686154, 0.09612142970736298,
        -0.05019554188579353, -0.00712019886112451],
       [0.011004476279509914, 0.0037363290518487783, 0.06898599713964727,
        0.09783996263076225, -0.008324136542912409],
       [0.0002600608894283072, 0.0007029296140874475,
        -0.0028374748790692653, 0.0030037276744491526,
        0.009190105262074236],
       [0.019504865390736998, -0.0729986176092966, -0.05861792972412624,
        -0.04164710074532667, -0.021473335179229933],
       [0.6211843273247705, -0.14379772284061648, 0.3284323141543874,
        0.5326401840557584, -0.24204522041622775],
       [-0.12262051441298684, 0.17799935256445415, -0.19239064287648583,
        -0.37276325724781895, 0.09286645244664224],
       [0.00884330651631477, -0.024538957579393735,
        0.0031269695223662443, 0.017399699740835755,
        -0.04145482498952329],
       [0.2431681327397587, 0.793447544402215, 0.6160934537812552,
        -0.1370748579009739, -0.025571966169774932],
       [-0.43998956774947473, 0.3318780691779483, -0.41690115331627675,
        0.3161678287889161, -0.010685873107643391],
       [0.01586454718058576, 0.07275137308761231, -0.006869903739187364,
        -0.007851463570964763, 0.03769196407806462]])

    df_pcs = pd.DataFrame(mat_pcs, columns=['PC1', 'PC2', 'PC3','PC4','PC5'], index=genes)

    return df_pcs

@pytest.fixture
def known_modes(genes):
    """
        known modes
    """
    mat_kns = np.array([[-0.14016657, -0.2187239 ,  0.01292835,  0.0528603 , -0.13300103,
        0.115241  ,  0.54783094,  0.03878433,  0.01508748, -0.15944897,
       -0.03015742, -0.03678603,  0.00142121, -0.1257714 ,  0.39561546,
       -0.19731769, -0.06508785,  0.29448432, -0.4854819 ,  0.17643273],
                   [-0.13137866,  0.03822342,  0.15881262,  0.01530638, -0.00873893,
       -0.28005314, -0.32090593,  0.47695143, -0.06338992,  0.15696721,
        0.21277554, -0.02011964,  0.19373257,  0.26277368, -0.28235581,
       -0.16583544, -0.04862901,  0.42497464, -0.16700115, -0.21875131]]).T

    df_kns = pd.DataFrame(mat_kns, columns=['EV1', 'EV2'], index=genes)

    return df_kns

@pytest.fixture
def expected_filter_result(genes):
    """
        expected filtered result
    """
    mat_result = np.array([[ 0.33187807, -0.41690115,  0.31616783, -0.01068587],
                            [-0.07299862, -0.05861793, -0.0416471 , -0.02147334],
                            [-0.03405309,  0.04666959,  0.01736164,  0.00555721],
                            [-0.10828103,  0.02153125,  0.04029951, -0.04690327],
                            [ 0.03914427,  0.30926583, -0.50488504,  0.1214284 ],
                            [ 0.07275137, -0.0068699 , -0.00785146,  0.03769196],
                            [-0.05426445,  0.39899045,  0.2535908 , -0.71596603],
                            [-0.06489108,  0.09612143, -0.05019554, -0.0071202 ],
                            [ 0.00070293, -0.00283747,  0.00300373,  0.00919011],
                            [ 0.00158758,  0.02828578, -0.02774233, -0.02366364],
                            [-0.02453896,  0.00312697,  0.0173997 , -0.04145482],
                            [ 0.00373633,  0.068986  ,  0.09783996, -0.00832414],
                            [ 0.39559131, -0.05559982,  0.34183307, -0.43065197],
                            [-0.12764441, -0.14143665, -0.02396365,  0.45733505],
                            [-0.03541021,  0.03249323,  0.0552054 ,  0.04305623],
                            [ 0.17799935, -0.19239064, -0.37276326,  0.09286645],
                            [ 0.00583199, -0.04226177, -0.04468407, -0.03976962],
                            [ 0.79344754,  0.61609345, -0.13707486, -0.02557197],
                            [-0.03240472,  0.00478096,  0.00308336, -0.0047098 ],
                            [-0.14379772,  0.32843231,  0.53264018, -0.24204522]])

    sorted_genes = ['6330527O06Rik', 'Arhgap36', 'Cd9', 'Cnp', 'Cnr1', 'Dmp1', 'Enpp2', 
                'Krt73', 'Lpar1', 'Mal', 'Mobp', 'Myl1', 'Npy', 'Penk', 'Plp1', 'Rgs4', 
                'Rgs5', 'Tac2', 'Trf', 'Vip']

    df_result = pd.DataFrame(mat_result, columns=['PC2', 'PC3', 'PC4', 'PC5'], index=sorted_genes)

    return df_result


def test_filter_known_modes(principal_components, known_modes, expected_filter_result):
    """
        test filter_known_modes function
    """

    df_pcs = principal_components
    df_kns = known_modes

    expected_result = expected_filter_result

    result = tc.filter_known_modes(df_pcs, df_kns)

    pd.testing.assert_frame_equal(result, expected_result)




