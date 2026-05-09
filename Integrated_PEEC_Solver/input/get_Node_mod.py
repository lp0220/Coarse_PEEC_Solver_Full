import math
import numpy as np


def read_node_file(node_file):
    """
    读取 Node 文件。

    文件格式：
    第一行：点的数量 N
    后面每行：x y z size
    """
    with open(node_file, "r") as f:
        lines = f.readlines()

    n_node = int(lines[0].strip())
    nodes = []

    for line in lines[1: n_node + 1]:
        parts = line.split()
        x = float(parts[0])
        y = float(parts[1])
        z = float(parts[2])
        size = float(parts[3])
        nodes.append([x, y, z, size])

    nodes = np.array(nodes, dtype=float)

    return nodes


def read_branch_file(branch_file):
    """
    读取 Branch 文件。

    文件格式：
    第一行：连接数量 M
    后面每行：node_i node_j

    注意：
    Branch 文件中的点编号通常是从 1 开始的，
    Python 数组索引是从 0 开始的，所以后面会减 1。
    """
    with open(branch_file, "r") as f:
        lines = f.readlines()

    n_branch = int(lines[0].strip())
    branches = []

    for line in lines[1: n_branch + 1]:
        parts = line.split()
        i = int(parts[0])
        j = int(parts[1])
        branches.append((i, j))

    return branches


def compute_new_dimension(nodes, branches, scale=1000.0):
    """
    计算每个点的新维度。

    对每个点：
    1. 找到所有连接点
    2. 坐标除以 scale
    3. 计算归一化后两点之间的欧氏距离 d
    4. 对 d 开根号
    5. 对所有相邻点结果求和
    6. 除以 sqrt(2)
    """
    n_node = nodes.shape[0]

    # 只取前三列坐标，并做归一化
    coords = nodes[:, 0:3] / scale

    # 每个点的新维度初始化为 0
    new_values = np.zeros(n_node, dtype=float)

    for i, j in branches:
        # 文件中编号是 1-based，转成 Python 的 0-based
        idx_i = i - 1
        idx_j = j - 1

        p_i = coords[idx_i]
        p_j = coords[idx_j]

        # 归一化坐标下的欧氏距离
        distance = np.linalg.norm(p_i - p_j)

        # 对距离开根号
        value = math.sqrt(distance)

        # 因为 i 和 j 是相互连接的，所以两个点都要加上这一项
        new_values[idx_i] += value
        new_values[idx_j] += value

    # 最后整体除以 sqrt(2)
    new_values = new_values / math.sqrt(1.0)

    return new_values


def write_output_file(output_file, nodes, new_values):
    """
    输出新的 Node 文件。

    格式：
    第一行：点数量
    后面每行：x y z size new_value
    """
    n_node = nodes.shape[0]

    with open(output_file, "w") as f:
        f.write(f"{n_node}\n")

        for i in range(n_node):
            x, y, z, size = nodes[i]
            new_value = new_values[i]

            f.write(
                f"{x:.16g} {y:.16g} {z:.16g} {size:.16g} {new_value:.16g}\n"
            )


def main():
    node_file = "Node.txt"
    branch_file = "Branch.txt"
    output_file = "Node_mod_sqrt1.txt"

    scale = 1000.0

    nodes = read_node_file(node_file)
    branches = read_branch_file(branch_file)

    new_values = compute_new_dimension(nodes, branches, scale)

    write_output_file(output_file, nodes, new_values)

    print("处理完成！")
    print(f"输入 Node 文件: {node_file}")
    print(f"输入 Branch 文件: {branch_file}")
    print(f"输出文件: {output_file}")
    print(f"scale = {scale}")

    # 打印前几个点的结果，方便检查
    print("\n前 10 个点的新维度结果：")
    for i in range(min(10, len(new_values))):
        print(f"Point {i + 1}: {new_values[i]:.12g}")


if __name__ == "__main__":
    main()