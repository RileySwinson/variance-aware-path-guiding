#pragma once
#if !defined(__DSCOMPARE_STRUCTURES_D_TREE_H_)
#define __DSCOMPARE_STRUCTURES_D_TREE_H_

#include <ds-compare/ds.h>
#include <Eigen/Core>

#include <array>
#include <atomic>
#include <chrono>
#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>

MTS_NAMESPACE_BEGIN

struct MTS_EXPORT_CORE DTreeRecord {
    Vector d;
    Float radiance, product;
    Float woPdf, bsdfPdf, dTreePdf;
    Float statisticalWeight;
    bool isDelta;
};

struct RandomGen {
    Float next1D() {
        return m_random->nextFloat();
    }

    Point2 next2D() {
        Float v1 = m_random->nextFloat();
        Float v2 = m_random->nextFloat();
        return Point2(v1, v2);
    }

private:
    ref<Random> m_random = new Random();
};

class MTS_EXPORT_CORE BlobWriter {
public:
    BlobWriter(const std::string& filename);

    template <typename Type>
    typename std::enable_if<std::is_standard_layout<Type>::value, BlobWriter&>::type operator<<(Type Element) {
        Write(&Element, 1);
        return *this;
    }

    // CAUTION: This function may break down on big-endian architectures.
    //          The ordering of bytes has to be reverted then.
    template <typename T>
    void Write(T* Src, size_t Size) {
        f.write(reinterpret_cast<const char*>(Src), Size * sizeof(T));
    }

private:
    std::ofstream f;
};

static void addToAtomicFloat(std::atomic<Float>& var, Float val) {
    auto current = var.load();
    while (!var.compare_exchange_weak(current, current + val));
}

class MTS_EXPORT_CORE AdamOptimizer {
public:
    AdamOptimizer(Float learningRate, int batchSize = 1, Float epsilon = 1e-08f, Float beta1 = 0.9f, Float beta2 = 0.999f);
    AdamOptimizer(const AdamOptimizer& arg);

    AdamOptimizer& operator=(const AdamOptimizer& arg);

    void append(Float gradient, Float statisticalWeight);
    void step(Float gradient);

    Float variable() const;

private:
    struct State {
        int iter = 0;
        Float firstMoment = 0;
        Float secondMoment = 0;
        Float variable = 0;

        Float batchAccumulation = 0;
        Float batchGradient = 0;
    } m_state;

    struct Hyperparameters {
        Float learningRate;
        int batchSize;
        Float epsilon;
        Float beta1;
        Float beta2;
    } m_hparams;
};

class MTS_EXPORT_CORE QuadTreeNode {
public:
    QuadTreeNode();
    QuadTreeNode(const QuadTreeNode& arg);

    QuadTreeNode& operator=(const QuadTreeNode& arg);

    void setSum(int index, Float val);
    void setSum(Float val);
    Float sum(int index) const;

    void setChild(int index, uint16_t val);
    uint16_t child(int index) const;

    void copyFrom(const QuadTreeNode& arg);
    int childIndex(Point2& p) const;
    
    // Evaluates the directional irradiance *sum density* (i.e. sum / area) at a given location p.
    // To obtain radiance, the sum density (result of this function) must be divided
    // by the total statistical weight of the estimates that were summed up.
    Float eval(Point2& p, const std::vector<QuadTreeNode>& nodes) const;
    
    Float pdf(Point2& p, const std::vector<QuadTreeNode>& nodes) const;
    int depthAt(Point2& p, const std::vector<QuadTreeNode>& nodes) const;
    Point2 sample(const std::vector<QuadTreeNode>& nodes) const;
    void record(Point2& p, Float irradiance, std::vector<QuadTreeNode>& nodes);
    void record(const Point2& origin, Float size, Point2 nodeOrigin, Float nodeSize, Float value, std::vector<QuadTreeNode>& nodes);
    Float computeOverlappingArea(const Point2& min1, const Point2& max1, const Point2& min2, const Point2& max2);
    bool isLeaf(int index) const;

    // Ensure that each quadtree node's sum of irradiance estimates equals that of all its children.
    void build(std::vector<QuadTreeNode>& nodes);

private:
    std::array<std::atomic<Float>, 4> m_sum;
    std::array<uint16_t, 4> m_children;
};

class MTS_EXPORT_CORE InternalDTree {
public:
    InternalDTree();

    const QuadTreeNode& node(size_t i) const;

    Float mean() const;

    void recordIrradiance(Point2 p, Float irradiance, Float statisticalWeight, DTreeParams::EDirectionalFilter directionalFilter);

    Float pdf(Point2 p) const;

    int depthAt(Point2 p) const;
    int depth() const;

    Point2 sample() const;

    size_t numNodes() const;

    Float statisticalWeight() const;
    void setStatisticalWeight(Float statisticalWeight);

    void reset(const InternalDTree& previousTree, int newMaxDepth, Float subdivisionThreshold);

    size_t approxMemoryFootprint() const;

    void build();

private:
    std::vector<QuadTreeNode> m_nodes;

    struct Atomic {
        Atomic() {
            sum.store(0, std::memory_order_relaxed);
            statisticalWeight.store(0, std::memory_order_relaxed);
        }

        Atomic(const Atomic& arg) {
            *this = arg;
        }

        Atomic& operator=(const Atomic& arg) {
            sum.store(arg.sum.load(std::memory_order_relaxed), std::memory_order_relaxed);
            statisticalWeight.store(arg.statisticalWeight.load(std::memory_order_relaxed), std::memory_order_relaxed);
            return *this;
        }

        std::atomic<Float> sum;
        std::atomic<Float> statisticalWeight;

    } m_atomic;

    int m_maxDepth;
};

struct MTS_EXPORT_CORE DirectionalTree : public DataStructure {
    static ref<RandomGen> random;

    ~DirectionalTree() { }

    void construct(DSArguments& init_data) override;

    void preprocess() override;

    void store(std::vector<Sample>& samples) override;

    void postprocess() override;

    Sample sample(Point2& pos) override;

    Float eval(Point2& pos) override;
    
    void wipe() override;

    DSType type() override;

    /* PPG funcs from here */

    void record(const DTreeRecord& rec, DTreeParams::EDirectionalFilter directionalFilter, DTreeParams::EBsdfSamplingFractionLoss bsdfSamplingFractionLoss);

    void build();

    void reset(int maxDepth, Float subdivisionThreshold);

    Vector sample() const;

    Float pdf(const Vector& dir) const;
    Float pdf(const Point2& uv) const;

    int depth() const;

    size_t numNodes() const;

    Float meanRadiance() const;

    Float statisticalWeight() const;

    Float statisticalWeightBuilding() const;

    void setStatisticalWeightBuilding(Float statisticalWeight);

    size_t approxMemoryFootprint() const;

    inline Float bsdfSamplingFraction(Float x) const {
        return 1 / (1 + std::exp(-x));
    }

    inline Float dBsdfSamplingFraction_dVariable(Float variable) const {
        Float fraction = bsdfSamplingFraction(variable);
        return fraction * (1 - fraction);
    }

    inline Float bsdfSamplingFraction() const {
        return bsdfSamplingFraction(bsdfSamplingFractionOptimizer.variable());
    }

    void optimizeBsdfSamplingFraction(const DTreeRecord& rec, Float ratioPower);

    void dump(BlobWriter& blob, const Point& p, const Vector& size) const;
private:
    InternalDTree building;
    InternalDTree sampling;
    std::vector<Sample> l_sample_storage;

    AdamOptimizer bsdfSamplingFractionOptimizer{0.01f};

    class SpinLock {
    public:
        SpinLock() {
            m_mutex.clear(std::memory_order_release);
        }

        SpinLock(const SpinLock& other) { m_mutex.clear(std::memory_order_release); }
        SpinLock& operator=(const SpinLock& other) { return *this; }

        void lock() {
            while (m_mutex.test_and_set(std::memory_order_acquire)) { }
        }

        void unlock() {
            m_mutex.clear(std::memory_order_release);
        }
    private:
        std::atomic_flag m_mutex;
    } m_lock;

    // Hyperparams (we're only interested in the directional ones)
    DTreeParams::EBsdfSamplingFractionLoss param_sampling_frac_loss;   // default = ENone
    DTreeParams::EDirectionalFilter        param_dir_filter;           // default = ENearest
    Float                                  param_d_tree_thresh;        // default = 0.01
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_STRUCTURES_D_TREE_H_ */