import mujoco
import scene_def as S


def test_mjcf_loads_and_has_geoms():
    path = S.generate_mjcf(S.SCENE)
    m = mujoco.MjModel.from_xml_path(path)
    names = {mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_GEOM, i) for i in range(m.ngeom)}
    assert "target" in names
    assert "obstacle0" in names
    assert "ws_ring" in names  # elliptical workspace outline
