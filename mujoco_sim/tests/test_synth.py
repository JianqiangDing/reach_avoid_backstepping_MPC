import textwrap
import scene_def as S
import synth


def test_scene_to_matlab_contains_sets_and_bound():
    src = synth.scene_to_matlab(S.SCENE)
    assert "safe_set_sym" in src and "target_set_sym" in src
    assert "a_max = 1.0" in src
    assert "y1" in src and "y2" in src
    assert "bound_min" in src and "bound_max" in src


def test_scene_hash_stable_and_sensitive():
    h1 = synth.scene_hash(S.SCENE)
    h2 = synth.scene_hash({**S.SCENE, "a_max": 2.0})
    assert h1 == synth.scene_hash(S.SCENE)
    assert h1 != h2


def test_load_bundle_from_fixture(tmp_path):
    mod = tmp_path / "fixture_ctrl.py"
    mod.write_text(textwrap.dedent('''
        from sympy import symbols
        x1, x2, x3, x4 = symbols("x1 x2 x3 x4")
        y1, y2 = symbols("y1 y2")
        u_opt = [x3, x4]
        certificate_opt = 1 - x1**2 - x2**2 - x3**2 - x4**2
        k1_opt = [-y1, -y2]
    '''))
    b = synth.load_bundle(str(mod))
    assert len(b["u_opt"]) == 2 and len(b["k1_opt"]) == 2
    assert b["certificate_opt"].free_symbols
