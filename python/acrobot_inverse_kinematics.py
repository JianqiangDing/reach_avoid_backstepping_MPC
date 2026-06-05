import matplotlib.pyplot as plt
import numpy as np


def kinematics(q1, q2, l1, l2):
    x = l1 * np.cos(q1) + l2 * np.cos(q1 + q2)
    y = l1 * np.sin(q1) + l2 * np.sin(q1 + q2)
    return (x, y)


def inverse_kinematics(y1, y2, l1, l2):
    """
    Compute inverse kinematics for a 2-link planar manipulator.

    Parameters:
        y1 (float): x-coordinate of end-effector
        y2 (float): y-coordinate of end-effector
        l1 (float): Length of first link
        l2 (float): Length of second link

    Returns:
        tuple: Two possible (q1, q2) angle solutions in radians
    """
    r2 = y1**2 + y2**2  # Squared distance to end-effector

    # Check if the target point is reachable
    # if r2 > (l1 + l2) ** 2 or r2 < (l1 - l2) ** 2:
    #     raise ValueError("Target point is outside the reachable workspace")

    # Compute q2 using the law of cosines
    cos_q2 = (r2 - l1**2 - l2**2) / (2 * l1 * l2)
    q2_1 = np.arccos(cos_q2)  # Elbow-down solution
    q2_2 = -np.arccos(cos_q2)  # Elbow-up solution

    # Compute q1 for both solutions
    phi = np.arctan2(y2, y1)
    k1_1 = l1 + l2 * np.cos(q2_1)
    k2_1 = l2 * np.sin(q2_1)
    q1_1 = phi - np.arctan2(k2_1, k1_1)

    k1_2 = l1 + l2 * np.cos(q2_2)
    k2_2 = l2 * np.sin(q2_2)
    q1_2 = phi - np.arctan2(k2_2, k1_2)

    return (q1_1, q2_1), (q1_2, q2_2)


if __name__ == "__main__":
    # Example usage
    y1, y2 = 1.0, -1.2  # Target position
    l1, l2 = 1.0, 1.0  # Link lengths

    (q1_sol1, q2_sol1), (q1_sol2, q2_sol2) = inverse_kinematics(y1, y2, l1, l2)

    (x_sol1, y_sol1) = kinematics(q1_sol1, q2_sol1, l1, l2)
    (x_sol2, y_sol2) = kinematics(q1_sol2, q2_sol2, l1, l2)

    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # First subplot
    axes[0].plot(
        [0, l1 * np.cos(q1_sol1)], [0, l1 * np.sin(q1_sol1)], "o-", lw=4, color="blue"
    )  # first link
    axes[0].plot(
        [l1 * np.cos(q1_sol1), x_sol1],
        [l1 * np.sin(q1_sol1), y_sol1],
        "o-",
        lw=4,
        color="red",
    )  # second link
    axes[0].set_xlim(-5, 5)
    axes[0].set_ylim(-5, 5)
    axes[0].axhline(0, color="black", linewidth=0.5)
    axes[0].axvline(0, color="black", linewidth=0.5)
    axes[0].grid()
    axes[0].set_title("Vector (x1, y1)")

    # Second subplot
    axes[1].plot(
        [0, l1 * np.cos(q1_sol2)], [0, l1 * np.sin(q1_sol2)], "o-", lw=4, color="blue"
    )  # first link
    axes[1].plot(
        [l1 * np.cos(q1_sol2), x_sol2],
        [l1 * np.sin(q1_sol2), y_sol2],
        "o-",
        lw=4,
        color="red",
    )  # second link
    axes[1].set_xlim(-5, 5)
    axes[1].set_ylim(-5, 5)
    axes[1].axhline(0, color="black", linewidth=0.5)
    axes[1].axvline(0, color="black", linewidth=0.5)
    axes[1].grid()
    axes[1].set_title("Vector (x2, y2)")

    # Show plot
    plt.show()
