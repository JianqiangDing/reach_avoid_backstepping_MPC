import numpy as np
import pytest
import scene_def as S


def test_safe_sign_inside_outside():
    sc = S.SCENE
    # a point in the workspace, clear of obstacles -> safe > 0
    assert S.safe_psi((0.54, 0.0), sc) > 0
    # inside obstacle 0 -> unsafe < 0
    ox, oy = sc["obstacles"][0]["center"]
    assert S.safe_psi((ox, oy), sc) < 0
    # outside the workspace ellipse -> unsafe < 0
    assert S.safe_psi((1.2, 0.0), sc) < 0


def test_target_sign():
    sc = S.SCENE
    tc = sc["target"]["center"]
    assert S.target_phi(tc, sc) < 0          # target center is inside target
    assert S.target_phi((0.54, 0.0), sc) > 0  # start side is not in target


def test_product_matches_truth_on_grid():
    sc = S.SCENE
    safe = S.safe_func(sc)
    cx, cy = sc["workspace"]["center"]; ax, ay = sc["workspace"]["semi"]
    gx, gy = np.meshgrid(np.linspace(cx - ax, cx + ax, 60),
                         np.linspace(cy - ay, cy + ay, 60))
    poly_safe = safe(gx, gy) >= 0
    in_ws = ((gx - cx) / ax) ** 2 + ((gy - cy) / ay) ** 2 <= 1
    out_obs = np.ones_like(in_ws, bool)
    for ob in sc["obstacles"]:
        ox, oy = ob["center"]
        out_obs &= (gx - ox) ** 2 + (gy - oy) ** 2 >= ob["radius"] ** 2
    truth = in_ws & out_obs
    # product sign agrees with the boolean safe region everywhere on the grid
    assert np.array_equal(poly_safe, truth)


def test_validate_rejects_overlap_and_exterior():
    assert S.validate_scene(S.SCENE) is True
    bad_overlap = {**S.SCENE, "obstacles": [
        {"center": (0.46, 0.0), "radius": 0.06},
        {"center": (0.49, 0.0), "radius": 0.06}]}
    with pytest.raises(ValueError):
        S.validate_scene(bad_overlap)
    bad_exterior = {**S.SCENE, "obstacles": [{"center": (0.46, 0.30), "radius": 0.05}]}
    with pytest.raises(ValueError):
        S.validate_scene(bad_exterior)
