import sympy as sp
import numpy as np
from scipy.integrate import odeint


def simulate(dyn, x0, dt, N=3):  # simulate one step from given x0
    if x0.ndim == 1:
        x0 = x0.reshape(1, -1)
    assert x0.ndim == 2

    def ode_fx(x, t):
        x = x.reshape(x0.shape)
        y = dyn(*x.T).squeeze(axis=1).T
        return y.flatten()

    x1 = odeint(ode_fx, x0.flatten(), np.linspace(0, dt, N))[-1].reshape(x0.shape)

    return x1


class BetterColor:
    @staticmethod
    def red0(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([213, 62, 79])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def red1(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([183, 35, 35])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def orange0(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([252, 141, 89])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def orange1(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([173, 108, 42])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def orange2(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([208, 186, 139])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def orange3(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([255, 140, 0])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def yellow0(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([254, 224, 139])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def yellow1(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([206, 190, 102])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def lightgreen(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([145, 207, 96])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def green0(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([26, 152, 80])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def green1(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([126, 197, 163])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def cyan0(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([50, 136, 189])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def cyan1(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([55, 146, 129])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def blue0(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([0.2, 0.3, 0.8]) * 255
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def blue1(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([114, 168, 206])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def purple0(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([84, 39, 136])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def purple1(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([167, 151, 210])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def gray0(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([178, 187, 205])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def gray1(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([201, 201, 201])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()

    @staticmethod
    def gray2(alpha=None):
        if alpha is not None:
            assert 0 <= alpha <= 1
        c = np.array([146, 149, 145])
        if alpha is None:
            return (c / 255).tolist()
        return np.concatenate([c / 255, [alpha]]).tolist()
