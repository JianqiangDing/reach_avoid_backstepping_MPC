import numpy as np
import scene_def as S
import verify_field as F


def test_field_reach_avoid_on_navigation_field():
    # a genuine reach-avoid navigation field: unit pull to target + short-range
    # obstacle repulsion (zero beyond an influence shell). field_reach_avoid should
    # report high reach and near-zero leave-safe for such a field.
    sc = S.SCENE
    tc = np.array(sc["target"]["center"])

    def k1(y):
        to_t = tc - y
        pull = to_t / (np.linalg.norm(to_t) + 1e-9)        # unit toward target
        rep = np.zeros(2)
        for ob in sc["obstacles"]:
            d = y - np.array(ob["center"]); dist = np.linalg.norm(d) + 1e-9
            influence = ob["radius"] + 0.06
            if dist < influence:
                rep += (influence - dist) / influence * d / dist * 3.0
        v = pull + rep
        return v / (np.linalg.norm(v) + 1e-9) * 0.3

    r = F.field_reach_avoid(sc, k1, n_starts=24, t_max=30.0, dt=0.02)
    assert r["success_frac"] > 0.8
    assert r["left_safe_frac"] < 0.2
