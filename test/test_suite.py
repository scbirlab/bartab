import numpy as np
from bartab.models.linear import WLSModel


def test_wls_uses_model_derived_weights():
    Y = np.array([
        [0.0, -0.2, -0.4, -0.6],
        [0.0, -0.2, -0.4, -0.6],
    ])
    x = np.array([0.0, 1.0, 2.0, 3.0])

    raw = np.array([
        [1000.0, 1000.0, 1000.0, 1000.0],  # reference
        [10.0, 20.0, 1000.0, 1000.0],      # variable uncertainty
    ])
    control_mask = np.array([True, False])
    groups = np.array(["a", "a", "b", "b"])

    model = WLSModel()
    weights = model.calculate_weights_matrix(
        Y,
        None,
        raw=raw,
        control_mask=control_mask,
        groups=groups,
    )

    assert weights.shape == Y.shape
    assert np.isfinite(weights).all()
    assert not np.allclose(weights, 1.0)
