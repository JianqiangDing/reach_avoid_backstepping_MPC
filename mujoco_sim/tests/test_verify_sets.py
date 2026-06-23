import scene_def as S
import verify_sets as V


def test_set_consistency_passes_for_default_scene():
    r = V.set_consistency(S.SCENE, n=120)
    assert r["agree"] is True
    assert r["mismatch_frac"] == 0.0
