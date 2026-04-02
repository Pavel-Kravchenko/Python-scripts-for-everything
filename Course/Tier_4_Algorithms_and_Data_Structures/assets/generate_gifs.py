import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch
import numpy as np
import os

# Create output directory
os.makedirs('gifs', exist_ok=True)

# Common settings
FIG_SIZE = (4, 3)  # ~400x300 pixels at 100 dpi
DPI = 100
FRAME_INTERVAL = 1200  # ms - slower for better comprehension


def create_bubble_sort_gif():
    """Create bubble sort animation showing swaps."""
    print("  Creating bubble-sort.gif...")
    
    arr = [64, 34, 25, 12, 22, 11, 90]
    n = len(arr)
    
    # Generate all states
    states = [(arr.copy(), -1, -1, "Initial array")]
    
    for i in range(n):
        for j in range(0, n - i - 1):
            states.append((arr.copy(), j, j + 1, f"Compare {arr[j]} and {arr[j+1]}"))
            if arr[j] > arr[j + 1]:
                arr[j], arr[j + 1] = arr[j + 1], arr[j]
                states.append((arr.copy(), j, j + 1, f"Swap! → {arr[j]}, {arr[j+1]}"))
    
    states.append((arr.copy(), -1, -1, "Sorted!"))
    
    fig, ax = plt.subplots(figsize=FIG_SIZE, dpi=DPI)
    
    def animate(frame):
        ax.clear()
        state, idx1, idx2, msg = states[frame % len(states)]
        
        colors = ['steelblue'] * len(state)
        if idx1 >= 0:
            colors[idx1] = 'red'
        if idx2 >= 0:
            colors[idx2] = 'orange'
        
        bars = ax.bar(range(len(state)), state, color=colors, edgecolor='black')
        
        for i, (bar, val) in enumerate(zip(bars, state)):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                   str(val), ha='center', va='bottom', fontsize=9, fontweight='bold')
        
        ax.set_xlim(-0.5, len(state) - 0.5)
        ax.set_ylim(0, max(state) + 15)
        ax.set_title(f'Bubble Sort\n{msg}', fontsize=10, fontweight='bold')
        ax.set_xticks([])
        ax.set_ylabel('Value')
        
    anim = animation.FuncAnimation(fig, animate, frames=len(states), interval=FRAME_INTERVAL)
    anim.save('gifs/bubble-sort.gif', writer='pillow')
    plt.close(fig)


def create_quicksort_gif():
    """Create quicksort animation with partition visualization."""
    print("  Creating quicksort.gif...")
    
    arr = [38, 27, 43, 3, 9, 82, 10]
    states = []
    
    def partition_visual(arr, low, high, states):
        pivot = arr[high]
        states.append((arr.copy(), high, -1, list(range(low, high)), f"Pivot = {pivot}"))
        
        i = low - 1
        for j in range(low, high):
            states.append((arr.copy(), high, j, list(range(low, high)), f"Compare {arr[j]} with pivot {pivot}"))
            if arr[j] <= pivot:
                i += 1
                arr[i], arr[j] = arr[j], arr[i]
                if i != j:
                    states.append((arr.copy(), high, -1, list(range(low, high)), f"Swap {arr[j]} and {arr[i]}"))
        
        arr[i + 1], arr[high] = arr[high], arr[i + 1]
        states.append((arr.copy(), i + 1, -1, [], f"Pivot placed at position {i+1}"))
        return i + 1
    
    def quicksort_visual(arr, low, high, states):
        if low < high:
            pi = partition_visual(arr, low, high, states)
            quicksort_visual(arr, low, pi - 1, states)
            quicksort_visual(arr, pi + 1, high, states)
    
    states.append((arr.copy(), -1, -1, [], "Initial array"))
    quicksort_visual(arr, 0, len(arr) - 1, states)
    states.append((arr.copy(), -1, -1, [], "Sorted!"))
    
    fig, ax = plt.subplots(figsize=FIG_SIZE, dpi=DPI)
    
    def animate(frame):
        ax.clear()
        state, pivot_idx, compare_idx, partition_range, msg = states[frame % len(states)]
        
        colors = ['steelblue'] * len(state)
        for i in partition_range:
            colors[i] = 'lightblue'
        if pivot_idx >= 0:
            colors[pivot_idx] = 'green'
        if compare_idx >= 0:
            colors[compare_idx] = 'red'
        
        bars = ax.bar(range(len(state)), state, color=colors, edgecolor='black')
        
        for i, (bar, val) in enumerate(zip(bars, state)):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                   str(val), ha='center', va='bottom', fontsize=9, fontweight='bold')
        
        ax.set_xlim(-0.5, len(state) - 0.5)
        ax.set_ylim(0, max(state) + 15)
        ax.set_title(f'QuickSort\n{msg}', fontsize=10, fontweight='bold')
        ax.set_xticks([])
        ax.set_ylabel('Value')
        
    anim = animation.FuncAnimation(fig, animate, frames=len(states), interval=FRAME_INTERVAL)
    anim.save('gifs/quicksort.gif', writer='pillow')
    plt.close(fig)


def create_merge_sort_gif():
    """Create merge sort animation showing divide and merge phases."""
    print("  Creating merge-sort.gif...")
    
    arr = [38, 27, 43, 3, 9, 82, 10]
    states = []
    
    def merge_sort_visual(arr, left, right, depth=0):
        if left < right:
            mid = (left + right) // 2
            states.append((arr.copy(), left, right, mid, f"Divide: [{left}:{right}] at mid={mid}"))
            
            merge_sort_visual(arr, left, mid, depth + 1)
            merge_sort_visual(arr, mid + 1, right, depth + 1)
            
            # Merge
            L = arr[left:mid + 1]
            R = arr[mid + 1:right + 1]
            i = j = 0
            k = left
            
            while i < len(L) and j < len(R):
                states.append((arr.copy(), left, right, mid, f"Merge: compare {L[i]} vs {R[j]}"))
                if L[i] <= R[j]:
                    arr[k] = L[i]
                    i += 1
                else:
                    arr[k] = R[j]
                    j += 1
                k += 1
            
            while i < len(L):
                arr[k] = L[i]
                i += 1
                k += 1
            
            while j < len(R):
                arr[k] = R[j]
                j += 1
                k += 1
            
            states.append((arr.copy(), left, right, mid, f"Merged [{left}:{right}]"))
    
    states.append((arr.copy(), 0, len(arr)-1, len(arr)//2, "Initial array"))
    merge_sort_visual(arr, 0, len(arr) - 1)
    states.append((arr.copy(), -1, -1, -1, "Sorted!"))
    
    fig, ax = plt.subplots(figsize=FIG_SIZE, dpi=DPI)
    
    def animate(frame):
        ax.clear()
        state, left, right, mid, msg = states[frame % len(states)]
        
        colors = ['steelblue'] * len(state)
        if left >= 0 and right >= 0:
            for i in range(left, right + 1):
                colors[i] = 'lightblue'
            if mid >= 0 and mid < len(state):
                colors[mid] = 'orange'
        
        bars = ax.bar(range(len(state)), state, color=colors, edgecolor='black')
        
        for i, (bar, val) in enumerate(zip(bars, state)):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                   str(val), ha='center', va='bottom', fontsize=9, fontweight='bold')
        
        ax.set_xlim(-0.5, len(state) - 0.5)
        ax.set_ylim(0, max(state) + 15)
        ax.set_title(f'Merge Sort\n{msg}', fontsize=10, fontweight='bold')
        ax.set_xticks([])
        ax.set_ylabel('Value')
        
    anim = animation.FuncAnimation(fig, animate, frames=len(states), interval=FRAME_INTERVAL)
    anim.save('gifs/merge-sort.gif', writer='pillow')
    plt.close(fig)


def create_binary_search_gif():
    """Create binary search animation narrowing down."""
    print("  Creating binary-search.gif...")
    
    arr = [2, 5, 8, 12, 16, 23, 38, 56, 72, 91]
    target = 23
    states = []
    
    left, right = 0, len(arr) - 1
    states.append((arr.copy(), left, right, -1, f"Find {target} in sorted array"))
    
    while left <= right:
        mid = (left + right) // 2
        states.append((arr.copy(), left, right, mid, f"Check mid={mid}: arr[{mid}]={arr[mid]}"))
        
        if arr[mid] == target:
            states.append((arr.copy(), mid, mid, mid, f"Found {target} at index {mid}!"))
            break
        elif arr[mid] < target:
            states.append((arr.copy(), left, right, mid, f"{arr[mid]} < {target}, go right"))
            left = mid + 1
        else:
            states.append((arr.copy(), left, right, mid, f"{arr[mid]} > {target}, go left"))
            right = mid - 1
    
    fig, ax = plt.subplots(figsize=FIG_SIZE, dpi=DPI)
    
    def animate(frame):
        ax.clear()
        state, left, right, mid, msg = states[frame % len(states)]
        
        colors = ['lightgray'] * len(state)
        for i in range(left, right + 1):
            colors[i] = 'steelblue'
        if mid >= 0:
            colors[mid] = 'red' if state[mid] != target else 'green'
        
        bars = ax.bar(range(len(state)), state, color=colors, edgecolor='black')
        
        for i, (bar, val) in enumerate(zip(bars, state)):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                   str(val), ha='center', va='bottom', fontsize=8, fontweight='bold')
        
        ax.set_xlim(-0.5, len(state) - 0.5)
        ax.set_ylim(0, max(state) + 15)
        ax.set_title(f'Binary Search (target={target})\n{msg}', fontsize=10, fontweight='bold')
        ax.set_xticks(range(len(state)))
        ax.set_xlabel('Index')
        ax.set_ylabel('Value')
        
    anim = animation.FuncAnimation(fig, animate, frames=len(states), interval=FRAME_INTERVAL)
    anim.save('gifs/binary-search.gif', writer='pillow')
    plt.close(fig)


def create_bst_insert_gif():
    """Create BST insertion animation."""
    print("  Creating bst-insert.gif...")
    
    class Node:
        def __init__(self, val):
            self.val = val
            self.left = None
            self.right = None
    
    def get_positions(node, x=0, y=0, dx=1.5, positions=None):
        if positions is None:
            positions = {}
        if node:
            positions[node.val] = (x, y)
            get_positions(node.left, x - dx, y - 1, dx * 0.6, positions)
            get_positions(node.right, x + dx, y - 1, dx * 0.6, positions)
        return positions
    
    def draw_tree(ax, node, positions, highlight=None, path=None):
        if path is None:
            path = []
        if node:
            x, y = positions[node.val]
            
            # Draw edges
            if node.left:
                lx, ly = positions[node.left.val]
                ax.plot([x, lx], [y, ly], 'k-', linewidth=2)
            if node.right:
                rx, ry = positions[node.right.val]
                ax.plot([x, rx], [y, ry], 'k-', linewidth=2)
            
            # Draw node
            color = 'green' if node.val == highlight else ('orange' if node.val in path else 'steelblue')
            circle = Circle((x, y), 0.3, color=color, ec='black', linewidth=2, zorder=5)
            ax.add_patch(circle)
            ax.text(x, y, str(node.val), ha='center', va='center', fontsize=10, 
                   fontweight='bold', color='white', zorder=6)
            
            draw_tree(ax, node.left, positions, highlight, path)
            draw_tree(ax, node.right, positions, highlight, path)
    
    values = [50, 30, 70, 20, 40, 60, 80]
    states = []
    
    root = None
    
    for val in values:
        path = []
        if root is None:
            root = Node(val)
            states.append((root, val, [], f"Insert {val} as root"))
        else:
            current = root
            while True:
                path.append(current.val)
                if val < current.val:
                    if current.left is None:
                        current.left = Node(val)
                        states.append((root, val, path.copy(), f"Insert {val} (left of {current.val})"))
                        break
                    current = current.left
                else:
                    if current.right is None:
                        current.right = Node(val)
                        states.append((root, val, path.copy(), f"Insert {val} (right of {current.val})"))
                        break
                    current = current.right
    
    states.append((root, None, [], "BST Complete!"))
    
    fig, ax = plt.subplots(figsize=FIG_SIZE, dpi=DPI)
    
    def animate(frame):
        ax.clear()
        tree, highlight, path, msg = states[frame % len(states)]
        
        positions = get_positions(tree)
        draw_tree(ax, tree, positions, highlight, path)
        
        ax.set_xlim(-3, 3)
        ax.set_ylim(-4, 1)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title(f'BST Insertion\n{msg}', fontsize=10, fontweight='bold')
        
    anim = animation.FuncAnimation(fig, animate, frames=len(states), interval=FRAME_INTERVAL)
    anim.save('gifs/bst-insert.gif', writer='pillow')
    plt.close(fig)


def create_avl_rotation_gif():
    """Create AVL tree rotation animation (single right rotation)."""
    print("  Creating avl-rotation.gif...")
    
    # We'll show a right rotation scenario
    # Before: y is root, x is left child, T2 is x's right subtree
    # After: x is root, y is right child, T2 is y's left subtree
    
    states = [
        ("unbalanced", "Unbalanced AVL Tree (Left-Left case)"),
        ("highlight_pivot", "Identify pivot node (y) for rotation"),
        ("rotating", "Performing Right Rotation..."),
        ("balanced", "Balanced! (Right Rotation Complete)")
    ]
    
    fig, ax = plt.subplots(figsize=FIG_SIZE, dpi=DPI)
    
    def draw_node(ax, x, y, val, color='steelblue'):
        circle = Circle((x, y), 0.35, color=color, ec='black', linewidth=2, zorder=5)
        ax.add_patch(circle)
        ax.text(x, y, str(val), ha='center', va='center', fontsize=11, 
               fontweight='bold', color='white', zorder=6)
    
    def draw_edge(ax, x1, y1, x2, y2):
        ax.plot([x1, x2], [y1, y2], 'k-', linewidth=2, zorder=1)
    
    def animate(frame):
        ax.clear()
        state, msg = states[frame % len(states)]
        
        if state == "unbalanced":
            # y at top, x below-left, z below-left of x
            draw_edge(ax, 0, 0, -1.2, -1)
            draw_edge(ax, 0, 0, 0.8, -1)
            draw_edge(ax, -1.2, -1, -1.8, -2)
            draw_edge(ax, -1.2, -1, -0.6, -2)
            
            draw_node(ax, 0, 0, 30, 'steelblue')
            draw_node(ax, -1.2, -1, 20, 'steelblue')
            draw_node(ax, 0.8, -1, 40, 'lightgray')
            draw_node(ax, -1.8, -2, 10, 'steelblue')
            draw_node(ax, -0.6, -2, 25, 'lightgray')
            
            ax.text(2, 0, "BF = -2", fontsize=9, color='red')
            
        elif state == "highlight_pivot":
            draw_edge(ax, 0, 0, -1.2, -1)
            draw_edge(ax, 0, 0, 0.8, -1)
            draw_edge(ax, -1.2, -1, -1.8, -2)
            draw_edge(ax, -1.2, -1, -0.6, -2)
            
            draw_node(ax, 0, 0, 30, 'red')  # pivot
            draw_node(ax, -1.2, -1, 20, 'orange')  # will become root
            draw_node(ax, 0.8, -1, 40, 'lightgray')
            draw_node(ax, -1.8, -2, 10, 'steelblue')
            draw_node(ax, -0.6, -2, 25, 'yellow')  # T2 - will move
            
        elif state == "rotating":
            # Intermediate state - showing movement
            draw_edge(ax, -0.6, 0, -1.5, -1)
            draw_edge(ax, -0.6, 0, 0.6, -1)
            draw_edge(ax, 0.6, -1, 0, -2)
            draw_edge(ax, 0.6, -1, 1.2, -2)
            
            draw_node(ax, -0.6, 0, 20, 'orange')
            draw_node(ax, -1.5, -1, 10, 'steelblue')
            draw_node(ax, 0.6, -1, 30, 'red')
            draw_node(ax, 0, -2, 25, 'yellow')
            draw_node(ax, 1.2, -2, 40, 'lightgray')
            
            ax.annotate('', xy=(0.3, 0.2), xytext=(-0.3, 0.5),
                       arrowprops=dict(arrowstyle='->', color='green', lw=2))
            
        elif state == "balanced":
            draw_edge(ax, 0, 0, -1.2, -1)
            draw_edge(ax, 0, 0, 1.2, -1)
            draw_edge(ax, 1.2, -1, 0.6, -2)
            draw_edge(ax, 1.2, -1, 1.8, -2)
            
            draw_node(ax, 0, 0, 20, 'green')
            draw_node(ax, -1.2, -1, 10, 'steelblue')
            draw_node(ax, 1.2, -1, 30, 'steelblue')
            draw_node(ax, 0.6, -2, 25, 'steelblue')
            draw_node(ax, 1.8, -2, 40, 'steelblue')
            
            ax.text(2, 0, "BF = 0", fontsize=9, color='green')
        
        ax.set_xlim(-3, 3)
        ax.set_ylim(-3, 1.5)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title(f'AVL Right Rotation\n{msg}', fontsize=10, fontweight='bold')
        
    anim = animation.FuncAnimation(fig, animate, frames=len(states), interval=FRAME_INTERVAL * 1.5)
    anim.save('gifs/avl-rotation.gif', writer='pillow')
    plt.close(fig)


def create_trie_insert_gif():
    """Create Trie construction animation word by word."""
    print("  Creating trie-insert.gif...")
    
    words = ["cat", "car", "card", "care"]
    
    class TrieNode:
        def __init__(self):
            self.children = {}
            self.is_end = False
    
    root = TrieNode()
    states = []
    
    def get_trie_positions(node, prefix="", x=0, y=0, dx=1.0, positions=None, edges=None):
        if positions is None:
            positions = {}
            edges = []
        
        positions[prefix if prefix else "root"] = (x, y, node.is_end)
        
        children = sorted(node.children.keys())
        n = len(children)
        
        for i, char in enumerate(children):
            child_x = x + (i - (n-1)/2) * dx
            child_y = y - 1
            child_prefix = prefix + char
            edges.append((prefix if prefix else "root", child_prefix, char))
            get_trie_positions(node.children[char], child_prefix, child_x, child_y, dx * 0.7, positions, edges)
        
        return positions, edges
    
    def copy_trie(node):
        new_node = TrieNode()
        new_node.is_end = node.is_end
        for char, child in node.children.items():
            new_node.children[char] = copy_trie(child)
        return new_node
    
    states.append((copy_trie(root), "", [], "Empty Trie"))
    
    for word in words:
        path = []
        current = root
        for i, char in enumerate(word):
            path.append(word[:i+1])
            if char not in current.children:
                current.children[char] = TrieNode()
                states.append((copy_trie(root), word, path.copy(), f"Insert '{word}': add '{char}'"))
            else:
                states.append((copy_trie(root), word, path.copy(), f"Insert '{word}': '{char}' exists"))
            current = current.children[char]
        current.is_end = True
        states.append((copy_trie(root), word, path.copy(), f"Mark end of '{word}'"))
    
    states.append((copy_trie(root), "", [], "Trie Complete!"))
    
    fig, ax = plt.subplots(figsize=FIG_SIZE, dpi=DPI)
    
    def animate(frame):
        ax.clear()
        trie, current_word, path, msg = states[frame % len(states)]
        
        positions, edges = get_trie_positions(trie)
        
        # Draw edges with labels
        for parent_key, child_key, char in edges:
            if parent_key in positions and child_key in positions:
                x1, y1, _ = positions[parent_key]
                x2, y2, _ = positions[child_key]
                ax.plot([x1, x2], [y1, y2], 'k-', linewidth=1.5)
                mx, my = (x1 + x2) / 2, (y1 + y2) / 2
                ax.text(mx + 0.1, my, char, fontsize=9, color='blue', fontweight='bold')
        
        # Draw nodes
        for key, (x, y, is_end) in positions.items():
            if key in path:
                color = 'orange'
            elif is_end:
                color = 'green'
            elif key == "root":
                color = 'steelblue'
            else:
                color = 'lightsteelblue'
            
            circle = Circle((x, y), 0.25, color=color, ec='black', linewidth=1.5, zorder=5)
            ax.add_patch(circle)
            
            if key == "root":
                ax.text(x, y, "∅", ha='center', va='center', fontsize=10, zorder=6)
        
        ax.set_xlim(-2.5, 2.5)
        ax.set_ylim(-4.5, 0.8)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title(f'Trie Construction\n{msg}', fontsize=10, fontweight='bold')
        
        # Legend
        ax.text(-2.3, 0.5, "● End of word", fontsize=8, color='green')
        
    anim = animation.FuncAnimation(fig, animate, frames=len(states), interval=FRAME_INTERVAL)
    anim.save('gifs/trie-insert.gif', writer='pillow')
    plt.close(fig)


if __name__ == "__main__":
    print("Generating algorithm GIFs...")
    print("=" * 40)
    
    create_bubble_sort_gif()
    create_quicksort_gif()
    create_merge_sort_gif()
    create_binary_search_gif()
    create_bst_insert_gif()
    create_avl_rotation_gif()
    create_trie_insert_gif()
    
    print("=" * 40)
    print("Done! GIFs saved to assets/gifs/")
