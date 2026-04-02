/**
 * Algorithms & Data Structures - Interactive Visualizations
 * Common JavaScript utilities and classes
 */

// ========== Utility Functions ==========

/**
 * Generate a random array of integers
 * @param {number} size - Array length
 * @param {number} min - Minimum value
 * @param {number} max - Maximum value
 * @returns {number[]} Random array
 */
function generateRandomArray(size, min = 5, max = 100) {
    const arr = [];
    for (let i = 0; i < size; i++) {
        arr.push(Math.floor(Math.random() * (max - min + 1)) + min);
    }
    return arr;
}

/**
 * Sleep for a specified duration
 * @param {number} ms - Milliseconds to sleep
 * @returns {Promise} Promise that resolves after ms
 */
function sleep(ms) {
    return new Promise(resolve => setTimeout(resolve, ms));
}

/**
 * Swap two elements in an array
 * @param {Array} arr - The array
 * @param {number} i - First index
 * @param {number} j - Second index
 */
function swap(arr, i, j) {
    [arr[i], arr[j]] = [arr[j], arr[i]];
}

/**
 * Parse custom array input string
 * @param {string} input - Comma-separated numbers
 * @returns {number[]} Parsed array
 */
function parseArrayInput(input) {
    return input.split(',')
        .map(s => parseInt(s.trim()))
        .filter(n => !isNaN(n) && n > 0 && n <= 100);
}

// ========== Visualization Base Class ==========

class AlgorithmVisualizer {
    constructor(containerId, options = {}) {
        this.container = document.getElementById(containerId);
        this.array = [];
        this.originalArray = [];
        this.isRunning = false;
        this.isPaused = false;
        this.isStepping = false;
        this.stepResolve = null;
        this.speed = options.speed || 500;
        this.comparisons = 0;
        this.swaps = 0;
        this.currentStep = 0;
        this.steps = [];
        
        this.onCompare = options.onCompare || (() => {});
        this.onSwap = options.onSwap || (() => {});
        this.onComplete = options.onComplete || (() => {});
    }

    /**
     * Initialize the array
     * @param {number[]} arr - Array to visualize
     */
    setArray(arr) {
        this.array = [...arr];
        this.originalArray = [...arr];
        this.reset();
    }

    /**
     * Reset statistics and state
     */
    reset() {
        this.comparisons = 0;
        this.swaps = 0;
        this.currentStep = 0;
        this.steps = [];
        this.isRunning = false;
        this.isPaused = false;
        this.updateStats();
    }

    /**
     * Render the array as bars
     */
    render(highlights = {}) {
        if (!this.container) return;
        
        const maxVal = Math.max(...this.array, 100);
        const containerHeight = this.container.clientHeight - 40;
        
        this.container.innerHTML = this.array.map((val, idx) => {
            const height = (val / maxVal) * containerHeight;
            let classes = ['bar'];
            
            if (highlights.comparing && highlights.comparing.includes(idx)) {
                classes.push('comparing');
            } else if (highlights.swapping && highlights.swapping.includes(idx)) {
                classes.push('swapping');
            } else if (highlights.sorted && highlights.sorted.includes(idx)) {
                classes.push('sorted');
            } else if (highlights.pivot === idx) {
                classes.push('pivot');
            } else if (highlights.current === idx) {
                classes.push('current');
            } else if (highlights.found === idx) {
                classes.push('found');
            } else if (highlights.eliminated && highlights.eliminated.includes(idx)) {
                classes.push('eliminated');
            }
            
            return `<div class="${classes.join(' ')}" style="height: ${height}px">${val}</div>`;
        }).join('');
    }

    /**
     * Update statistics display
     */
    updateStats() {
        const compEl = document.getElementById('comparisons');
        const swapEl = document.getElementById('swaps');
        const stepEl = document.getElementById('current-step');
        
        if (compEl) compEl.textContent = this.comparisons;
        if (swapEl) swapEl.textContent = this.swaps;
        if (stepEl) stepEl.textContent = this.currentStep;
    }

    /**
     * Add a step to the history
     */
    addStep(description) {
        this.currentStep++;
        this.steps.push({
            number: this.currentStep,
            description,
            array: [...this.array]
        });
        this.updateStepsList();
        this.updateStats();
    }

    /**
     * Update the steps list display
     */
    updateStepsList() {
        const list = document.getElementById('steps-list');
        if (!list) return;
        
        const step = this.steps[this.steps.length - 1];
        const stepEl = document.createElement('div');
        stepEl.className = 'step current';
        stepEl.innerHTML = `<span class="step-number">${step.number}.</span> ${step.description}`;
        
        // Remove current class from previous step
        const prev = list.querySelector('.step.current');
        if (prev) prev.classList.remove('current');
        
        list.appendChild(stepEl);
        list.scrollTop = list.scrollHeight;
    }

    /**
     * Clear steps list
     */
    clearSteps() {
        const list = document.getElementById('steps-list');
        if (list) list.innerHTML = '';
    }

    /**
     * Wait for animation delay or step button
     */
    async wait() {
        if (this.isStepping) {
            await new Promise(resolve => {
                this.stepResolve = resolve;
            });
        } else {
            await sleep(this.speed);
            while (this.isPaused && !this.isStepping) {
                await sleep(100);
            }
        }
    }

    /**
     * Handle step button click
     */
    step() {
        if (this.stepResolve) {
            this.stepResolve();
            this.stepResolve = null;
        }
    }

    /**
     * Pause execution
     */
    pause() {
        this.isPaused = true;
    }

    /**
     * Resume execution
     */
    resume() {
        this.isPaused = false;
    }

    /**
     * Stop execution
     */
    stop() {
        this.isRunning = false;
        this.isPaused = false;
        this.isStepping = false;
        if (this.stepResolve) {
            this.stepResolve();
        }
    }

    /**
     * Set animation speed
     */
    setSpeed(speed) {
        this.speed = speed;
    }

    /**
     * Highlight code line
     */
    highlightCode(lineNumber) {
        const lines = document.querySelectorAll('.code-line');
        lines.forEach((line, idx) => {
            line.classList.toggle('highlight', idx + 1 === lineNumber);
        });
    }
}

// ========== Sorting Algorithms ==========

class SortingVisualizer extends AlgorithmVisualizer {
    constructor(containerId, options = {}) {
        super(containerId, options);
        this.sortedIndices = new Set();
    }

    reset() {
        super.reset();
        this.sortedIndices = new Set();
    }

    /**
     * Bubble Sort implementation
     */
    async bubbleSort() {
        this.isRunning = true;
        this.clearSteps();
        const n = this.array.length;
        
        for (let i = 0; i < n - 1 && this.isRunning; i++) {
            let swapped = false;
            
            for (let j = 0; j < n - i - 1 && this.isRunning; j++) {
                // Compare
                this.comparisons++;
                this.render({ comparing: [j, j + 1], sorted: [...this.sortedIndices] });
                this.addStep(`Compare arr[${j}]=${this.array[j]} with arr[${j+1}]=${this.array[j+1]}`);
                this.highlightCode(5);
                await this.wait();
                
                if (this.array[j] > this.array[j + 1]) {
                    // Swap
                    this.swaps++;
                    swap(this.array, j, j + 1);
                    this.render({ swapping: [j, j + 1], sorted: [...this.sortedIndices] });
                    this.addStep(`Swap: ${this.array[j+1]} ↔ ${this.array[j]}`);
                    this.highlightCode(6);
                    await this.wait();
                    swapped = true;
                }
            }
            
            // Mark last element as sorted
            this.sortedIndices.add(n - i - 1);
            this.render({ sorted: [...this.sortedIndices] });
            
            if (!swapped) {
                this.addStep('No swaps needed - array is sorted!');
                break;
            }
        }
        
        // Mark all as sorted
        for (let i = 0; i < n; i++) this.sortedIndices.add(i);
        this.render({ sorted: [...this.sortedIndices] });
        this.highlightCode(0);
        this.isRunning = false;
        this.onComplete();
    }

    /**
     * Merge Sort implementation
     */
    async mergeSort() {
        this.isRunning = true;
        this.clearSteps();
        await this._mergeSortHelper(0, this.array.length - 1);
        
        // Mark all as sorted
        for (let i = 0; i < this.array.length; i++) this.sortedIndices.add(i);
        this.render({ sorted: [...this.sortedIndices] });
        this.isRunning = false;
        this.onComplete();
    }

    async _mergeSortHelper(left, right) {
        if (left >= right || !this.isRunning) return;
        
        const mid = Math.floor((left + right) / 2);
        this.addStep(`Divide: [${left}..${mid}] and [${mid+1}..${right}]`);
        this.highlightCode(3);
        
        await this._mergeSortHelper(left, mid);
        await this._mergeSortHelper(mid + 1, right);
        await this._merge(left, mid, right);
    }

    async _merge(left, mid, right) {
        if (!this.isRunning) return;
        
        const leftArr = this.array.slice(left, mid + 1);
        const rightArr = this.array.slice(mid + 1, right + 1);
        
        this.addStep(`Merge: [${leftArr.join(',')}] and [${rightArr.join(',')}]`);
        this.highlightCode(7);
        
        let i = 0, j = 0, k = left;
        
        while (i < leftArr.length && j < rightArr.length && this.isRunning) {
            this.comparisons++;
            this.render({ comparing: [left + i, mid + 1 + j] });
            await this.wait();
            
            if (leftArr[i] <= rightArr[j]) {
                this.array[k] = leftArr[i];
                i++;
            } else {
                this.array[k] = rightArr[j];
                j++;
            }
            this.swaps++;
            this.render({ swapping: [k] });
            await this.wait();
            k++;
        }
        
        while (i < leftArr.length && this.isRunning) {
            this.array[k] = leftArr[i];
            i++;
            k++;
        }
        
        while (j < rightArr.length && this.isRunning) {
            this.array[k] = rightArr[j];
            j++;
            k++;
        }
        
        this.render();
    }

    /**
     * QuickSort implementation
     */
    async quickSort() {
        this.isRunning = true;
        this.clearSteps();
        await this._quickSortHelper(0, this.array.length - 1);
        
        // Mark all as sorted
        for (let i = 0; i < this.array.length; i++) this.sortedIndices.add(i);
        this.render({ sorted: [...this.sortedIndices] });
        this.isRunning = false;
        this.onComplete();
    }

    async _quickSortHelper(low, high) {
        if (low >= high || !this.isRunning) return;
        
        const pivotIdx = await this._partition(low, high);
        this.sortedIndices.add(pivotIdx);
        
        await this._quickSortHelper(low, pivotIdx - 1);
        await this._quickSortHelper(pivotIdx + 1, high);
    }

    async _partition(low, high) {
        const pivot = this.array[high];
        this.addStep(`Pivot selected: ${pivot} at index ${high}`);
        this.highlightCode(4);
        this.render({ pivot: high });
        await this.wait();
        
        let i = low - 1;
        
        for (let j = low; j < high && this.isRunning; j++) {
            this.comparisons++;
            this.render({ comparing: [j], pivot: high, sorted: [...this.sortedIndices] });
            this.addStep(`Compare arr[${j}]=${this.array[j]} with pivot ${pivot}`);
            this.highlightCode(7);
            await this.wait();
            
            if (this.array[j] < pivot) {
                i++;
                if (i !== j) {
                    this.swaps++;
                    swap(this.array, i, j);
                    this.render({ swapping: [i, j], pivot: high, sorted: [...this.sortedIndices] });
                    this.addStep(`Swap arr[${i}] ↔ arr[${j}]`);
                    this.highlightCode(9);
                    await this.wait();
                }
            }
        }
        
        // Place pivot in correct position
        this.swaps++;
        swap(this.array, i + 1, high);
        this.render({ swapping: [i + 1, high], sorted: [...this.sortedIndices] });
        this.addStep(`Place pivot at index ${i + 1}`);
        this.highlightCode(12);
        await this.wait();
        
        return i + 1;
    }

    /**
     * Selection Sort implementation
     */
    async selectionSort() {
        this.isRunning = true;
        this.clearSteps();
        const n = this.array.length;
        
        for (let i = 0; i < n - 1 && this.isRunning; i++) {
            let minIdx = i;
            this.addStep(`Find minimum in unsorted portion [${i}..${n-1}]`);
            this.highlightCode(3);
            
            for (let j = i + 1; j < n && this.isRunning; j++) {
                this.comparisons++;
                this.render({ comparing: [minIdx, j], current: i, sorted: [...this.sortedIndices] });
                await this.wait();
                
                if (this.array[j] < this.array[minIdx]) {
                    minIdx = j;
                    this.addStep(`New minimum found: ${this.array[minIdx]} at index ${minIdx}`);
                }
            }
            
            if (minIdx !== i) {
                this.swaps++;
                swap(this.array, i, minIdx);
                this.render({ swapping: [i, minIdx], sorted: [...this.sortedIndices] });
                this.addStep(`Swap arr[${i}] ↔ arr[${minIdx}]`);
                this.highlightCode(8);
                await this.wait();
            }
            
            this.sortedIndices.add(i);
        }
        
        this.sortedIndices.add(n - 1);
        this.render({ sorted: [...this.sortedIndices] });
        this.isRunning = false;
        this.onComplete();
    }

    /**
     * Insertion Sort implementation
     */
    async insertionSort() {
        this.isRunning = true;
        this.clearSteps();
        const n = this.array.length;
        this.sortedIndices.add(0);
        
        for (let i = 1; i < n && this.isRunning; i++) {
            const key = this.array[i];
            let j = i - 1;
            
            this.addStep(`Insert ${key} into sorted portion [0..${i-1}]`);
            this.highlightCode(3);
            this.render({ current: i, sorted: [...this.sortedIndices] });
            await this.wait();
            
            while (j >= 0 && this.array[j] > key && this.isRunning) {
                this.comparisons++;
                this.render({ comparing: [j, j + 1], sorted: [...this.sortedIndices] });
                await this.wait();
                
                this.array[j + 1] = this.array[j];
                this.swaps++;
                this.render({ swapping: [j, j + 1] });
                this.addStep(`Shift ${this.array[j]} to position ${j + 1}`);
                this.highlightCode(6);
                await this.wait();
                j--;
            }
            
            this.array[j + 1] = key;
            this.sortedIndices.add(i);
            this.render({ sorted: [...this.sortedIndices] });
        }
        
        this.isRunning = false;
        this.onComplete();
    }
}

// ========== Searching Algorithms ==========

class SearchingVisualizer extends AlgorithmVisualizer {
    constructor(containerId, options = {}) {
        super(containerId, options);
        this.target = null;
        this.foundIndex = -1;
    }

    setTarget(target) {
        this.target = target;
        this.foundIndex = -1;
    }

    /**
     * Linear Search implementation
     */
    async linearSearch() {
        this.isRunning = true;
        this.clearSteps();
        this.foundIndex = -1;
        
        for (let i = 0; i < this.array.length && this.isRunning; i++) {
            this.comparisons++;
            this.render({ current: i });
            this.addStep(`Check arr[${i}] = ${this.array[i]} against target ${this.target}`);
            this.highlightCode(3);
            await this.wait();
            
            if (this.array[i] === this.target) {
                this.foundIndex = i;
                this.render({ found: i });
                this.addStep(`Found! Target ${this.target} at index ${i}`);
                this.highlightCode(4);
                this.isRunning = false;
                this.onComplete();
                return i;
            }
        }
        
        this.addStep(`Target ${this.target} not found in array`);
        this.isRunning = false;
        this.onComplete();
        return -1;
    }

    /**
     * Binary Search implementation
     */
    async binarySearch() {
        this.isRunning = true;
        this.clearSteps();
        this.foundIndex = -1;
        
        let left = 0;
        let right = this.array.length - 1;
        const eliminated = new Set();
        
        while (left <= right && this.isRunning) {
            const mid = Math.floor((left + right) / 2);
            
            this.comparisons++;
            this.render({ current: mid, eliminated: [...eliminated] });
            this.addStep(`Check mid=${mid}: arr[${mid}] = ${this.array[mid]}`);
            this.highlightCode(4);
            await this.wait();
            
            if (this.array[mid] === this.target) {
                this.foundIndex = mid;
                this.render({ found: mid });
                this.addStep(`Found! Target ${this.target} at index ${mid}`);
                this.highlightCode(5);
                this.isRunning = false;
                this.onComplete();
                return mid;
            } else if (this.array[mid] < this.target) {
                // Eliminate left half
                for (let i = left; i <= mid; i++) eliminated.add(i);
                this.addStep(`${this.array[mid]} < ${this.target}, search right half`);
                this.highlightCode(7);
                left = mid + 1;
            } else {
                // Eliminate right half
                for (let i = mid; i <= right; i++) eliminated.add(i);
                this.addStep(`${this.array[mid]} > ${this.target}, search left half`);
                this.highlightCode(9);
                right = mid - 1;
            }
            
            this.render({ eliminated: [...eliminated] });
            await this.wait();
        }
        
        this.addStep(`Target ${this.target} not found in array`);
        this.isRunning = false;
        this.onComplete();
        return -1;
    }
}

// ========== Tree Visualization ==========

class TreeNode {
    constructor(value) {
        this.value = value;
        this.left = null;
        this.right = null;
        this.height = 1; // For AVL
        this.x = 0;
        this.y = 0;
    }
}

class BSTVisualizer {
    constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        this.ctx = this.canvas ? this.canvas.getContext('2d') : null;
        this.root = null;
        this.nodeRadius = 25;
        this.levelHeight = 70;
        this.speed = 500;
        this.isRunning = false;
        this.highlightedNodes = new Set();
        this.foundNode = null;
        this.comparisons = 0;
        this.steps = [];
    }

    reset() {
        this.root = null;
        this.highlightedNodes = new Set();
        this.foundNode = null;
        this.comparisons = 0;
        this.steps = [];
        this.render();
    }

    /**
     * Calculate positions for all nodes
     */
    calculatePositions() {
        if (!this.root || !this.canvas) return;
        
        const width = this.canvas.width;
        this._calculatePositions(this.root, width / 2, 50, width / 4);
    }

    _calculatePositions(node, x, y, offset) {
        if (!node) return;
        
        node.x = x;
        node.y = y;
        
        if (node.left) {
            this._calculatePositions(node.left, x - offset, y + this.levelHeight, offset / 2);
        }
        if (node.right) {
            this._calculatePositions(node.right, x + offset, y + this.levelHeight, offset / 2);
        }
    }

    /**
     * Render the tree
     */
    render() {
        if (!this.ctx || !this.canvas) return;
        
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
        this.calculatePositions();
        
        if (this.root) {
            this._drawEdges(this.root);
            this._drawNodes(this.root);
        }
    }

    _drawEdges(node) {
        if (!node) return;
        
        this.ctx.strokeStyle = '#7f8c8d';
        this.ctx.lineWidth = 2;
        
        if (node.left) {
            this.ctx.beginPath();
            this.ctx.moveTo(node.x, node.y);
            this.ctx.lineTo(node.left.x, node.left.y);
            this.ctx.stroke();
            this._drawEdges(node.left);
        }
        
        if (node.right) {
            this.ctx.beginPath();
            this.ctx.moveTo(node.x, node.y);
            this.ctx.lineTo(node.right.x, node.right.y);
            this.ctx.stroke();
            this._drawEdges(node.right);
        }
    }

    _drawNodes(node) {
        if (!node) return;
        
        // Draw node circle
        this.ctx.beginPath();
        this.ctx.arc(node.x, node.y, this.nodeRadius, 0, Math.PI * 2);
        
        if (this.foundNode === node) {
            this.ctx.fillStyle = '#2ecc71';
            this.ctx.strokeStyle = '#27ae60';
        } else if (this.highlightedNodes.has(node)) {
            this.ctx.fillStyle = '#3498db';
            this.ctx.strokeStyle = '#2980b9';
        } else {
            this.ctx.fillStyle = '#3498db';
            this.ctx.strokeStyle = '#2980b9';
        }
        
        this.ctx.fill();
        this.ctx.lineWidth = 3;
        this.ctx.stroke();
        
        // Draw value
        this.ctx.fillStyle = 'white';
        this.ctx.font = 'bold 16px Arial';
        this.ctx.textAlign = 'center';
        this.ctx.textBaseline = 'middle';
        this.ctx.fillText(node.value.toString(), node.x, node.y);
        
        this._drawNodes(node.left);
        this._drawNodes(node.right);
    }

    /**
     * Insert a value into the BST
     */
    async insert(value) {
        this.isRunning = true;
        this.highlightedNodes.clear();
        this.foundNode = null;
        
        if (!this.root) {
            this.root = new TreeNode(value);
            this.addStep(`Insert ${value} as root`);
            this.render();
            this.isRunning = false;
            return;
        }
        
        await this._insertHelper(this.root, value);
        this.highlightedNodes.clear();
        this.render();
        this.isRunning = false;
    }

    async _insertHelper(node, value) {
        this.comparisons++;
        this.highlightedNodes.add(node);
        this.render();
        this.addStep(`Compare ${value} with ${node.value}`);
        await sleep(this.speed);
        
        if (value < node.value) {
            if (node.left) {
                await this._insertHelper(node.left, value);
            } else {
                node.left = new TreeNode(value);
                this.addStep(`Insert ${value} as left child of ${node.value}`);
            }
        } else {
            if (node.right) {
                await this._insertHelper(node.right, value);
            } else {
                node.right = new TreeNode(value);
                this.addStep(`Insert ${value} as right child of ${node.value}`);
            }
        }
    }

    /**
     * Search for a value in the BST
     */
    async search(value) {
        this.isRunning = true;
        this.highlightedNodes.clear();
        this.foundNode = null;
        this.comparisons = 0;
        
        const result = await this._searchHelper(this.root, value);
        this.isRunning = false;
        return result;
    }

    async _searchHelper(node, value) {
        if (!node) {
            this.addStep(`Value ${value} not found`);
            return null;
        }
        
        this.comparisons++;
        this.highlightedNodes.add(node);
        this.render();
        this.addStep(`Compare ${value} with ${node.value}`);
        await sleep(this.speed);
        
        if (value === node.value) {
            this.foundNode = node;
            this.render();
            this.addStep(`Found ${value}!`);
            return node;
        } else if (value < node.value) {
            this.addStep(`${value} < ${node.value}, go left`);
            return await this._searchHelper(node.left, value);
        } else {
            this.addStep(`${value} > ${node.value}, go right`);
            return await this._searchHelper(node.right, value);
        }
    }

    /**
     * Delete a value from the BST
     */
    async delete(value) {
        this.isRunning = true;
        this.highlightedNodes.clear();
        this.foundNode = null;
        
        this.root = await this._deleteHelper(this.root, value);
        this.highlightedNodes.clear();
        this.render();
        this.isRunning = false;
    }

    async _deleteHelper(node, value) {
        if (!node) {
            this.addStep(`Value ${value} not found`);
            return null;
        }
        
        this.highlightedNodes.add(node);
        this.render();
        await sleep(this.speed);
        
        if (value < node.value) {
            node.left = await this._deleteHelper(node.left, value);
        } else if (value > node.value) {
            node.right = await this._deleteHelper(node.right, value);
        } else {
            // Node found
            this.addStep(`Found ${value}, deleting...`);
            
            if (!node.left) {
                return node.right;
            } else if (!node.right) {
                return node.left;
            }
            
            // Node with two children
            let successor = node.right;
            while (successor.left) {
                successor = successor.left;
            }
            
            this.addStep(`Replace with successor ${successor.value}`);
            node.value = successor.value;
            node.right = await this._deleteHelper(node.right, successor.value);
        }
        
        return node;
    }

    addStep(description) {
        this.steps.push(description);
        this.updateStepsList();
    }

    updateStepsList() {
        const list = document.getElementById('steps-list');
        if (!list) return;
        
        const step = this.steps[this.steps.length - 1];
        const stepEl = document.createElement('div');
        stepEl.className = 'step current';
        stepEl.innerHTML = `<span class="step-number">${this.steps.length}.</span> ${step}`;
        
        const prev = list.querySelector('.step.current');
        if (prev) prev.classList.remove('current');
        
        list.appendChild(stepEl);
        list.scrollTop = list.scrollHeight;
    }

    clearSteps() {
        this.steps = [];
        const list = document.getElementById('steps-list');
        if (list) list.innerHTML = '';
    }

    setSpeed(speed) {
        this.speed = speed;
    }
}

// ========== AVL Tree ==========

class AVLVisualizer extends BSTVisualizer {
    constructor(canvasId) {
        super(canvasId);
        this.rotationSteps = [];
    }

    getHeight(node) {
        return node ? node.height : 0;
    }

    getBalance(node) {
        return node ? this.getHeight(node.left) - this.getHeight(node.right) : 0;
    }

    updateHeight(node) {
        if (node) {
            node.height = 1 + Math.max(this.getHeight(node.left), this.getHeight(node.right));
        }
    }

    async rightRotate(y) {
        this.addStep(`Right Rotate at ${y.value}`);
        await sleep(this.speed);
        
        const x = y.left;
        const T2 = x.right;
        
        x.right = y;
        y.left = T2;
        
        this.updateHeight(y);
        this.updateHeight(x);
        
        this.render();
        await sleep(this.speed);
        
        return x;
    }

    async leftRotate(x) {
        this.addStep(`Left Rotate at ${x.value}`);
        await sleep(this.speed);
        
        const y = x.right;
        const T2 = y.left;
        
        y.left = x;
        x.right = T2;
        
        this.updateHeight(x);
        this.updateHeight(y);
        
        this.render();
        await sleep(this.speed);
        
        return y;
    }

    async insert(value) {
        this.isRunning = true;
        this.highlightedNodes.clear();
        this.clearSteps();
        
        this.root = await this._avlInsert(this.root, value);
        
        this.highlightedNodes.clear();
        this.render();
        this.isRunning = false;
    }

    async _avlInsert(node, value) {
        if (!node) {
            this.addStep(`Insert ${value}`);
            return new TreeNode(value);
        }
        
        this.highlightedNodes.add(node);
        this.render();
        this.addStep(`Compare ${value} with ${node.value}`);
        await sleep(this.speed);
        
        if (value < node.value) {
            node.left = await this._avlInsert(node.left, value);
        } else if (value > node.value) {
            node.right = await this._avlInsert(node.right, value);
        } else {
            return node; // Duplicate
        }
        
        this.updateHeight(node);
        const balance = this.getBalance(node);
        
        this.addStep(`Balance factor at ${node.value}: ${balance}`);
        this.render();
        await sleep(this.speed / 2);
        
        // Left Left
        if (balance > 1 && value < node.left.value) {
            this.addStep(`LL Case detected at ${node.value}`);
            return await this.rightRotate(node);
        }
        
        // Right Right
        if (balance < -1 && value > node.right.value) {
            this.addStep(`RR Case detected at ${node.value}`);
            return await this.leftRotate(node);
        }
        
        // Left Right
        if (balance > 1 && value > node.left.value) {
            this.addStep(`LR Case detected at ${node.value}`);
            node.left = await this.leftRotate(node.left);
            return await this.rightRotate(node);
        }
        
        // Right Left
        if (balance < -1 && value < node.right.value) {
            this.addStep(`RL Case detected at ${node.value}`);
            node.right = await this.rightRotate(node.right);
            return await this.leftRotate(node);
        }
        
        return node;
    }

    _drawNodes(node) {
        if (!node) return;
        
        // Draw node circle
        this.ctx.beginPath();
        this.ctx.arc(node.x, node.y, this.nodeRadius, 0, Math.PI * 2);
        
        if (this.foundNode === node) {
            this.ctx.fillStyle = '#2ecc71';
        } else if (this.highlightedNodes.has(node)) {
            this.ctx.fillStyle = '#e74c3c';
        } else {
            this.ctx.fillStyle = '#3498db';
        }
        
        this.ctx.fill();
        this.ctx.strokeStyle = '#2c3e50';
        this.ctx.lineWidth = 3;
        this.ctx.stroke();
        
        // Draw value
        this.ctx.fillStyle = 'white';
        this.ctx.font = 'bold 14px Arial';
        this.ctx.textAlign = 'center';
        this.ctx.textBaseline = 'middle';
        this.ctx.fillText(node.value.toString(), node.x, node.y - 5);
        
        // Draw height/balance
        this.ctx.font = '10px Arial';
        this.ctx.fillText(`h:${node.height}`, node.x, node.y + 10);
        
        this._drawNodes(node.left);
        this._drawNodes(node.right);
    }
}

// ========== Export for use in HTML ==========
if (typeof window !== 'undefined') {
    window.AlgorithmVisualizer = AlgorithmVisualizer;
    window.SortingVisualizer = SortingVisualizer;
    window.SearchingVisualizer = SearchingVisualizer;
    window.BSTVisualizer = BSTVisualizer;
    window.AVLVisualizer = AVLVisualizer;
    window.generateRandomArray = generateRandomArray;
    window.parseArrayInput = parseArrayInput;
    window.sleep = sleep;
}
