#include <mitsuba/ds/structures/d_tree.h>

MTS_NAMESPACE_BEGIN

/* ========== */
/* BlobWriter */
/* ========== */

BlobWriter::BlobWriter(const std::string& filename) : f(filename, std::ios::out | std::ios::binary) { }

/* ============= */
/* AdamOptimizer */
/* ============= */

AdamOptimizer::AdamOptimizer(Float learningRate, int batchSize, Float epsilon, Float beta1, Float beta2) {
	m_hparams = { learningRate, batchSize, epsilon, beta1, beta2 };
}

AdamOptimizer::AdamOptimizer(const AdamOptimizer& arg) {
    *this = arg;
}

AdamOptimizer& AdamOptimizer::operator=(const AdamOptimizer& arg) {
    m_state = arg.m_state;
    m_hparams = arg.m_hparams;
    return *this;
}

void AdamOptimizer::append(Float gradient, Float statisticalWeight) {
    m_state.batchGradient += gradient * statisticalWeight;
    m_state.batchAccumulation += statisticalWeight;

    if (m_state.batchAccumulation > m_hparams.batchSize) {
        step(m_state.batchGradient / m_state.batchAccumulation);

        m_state.batchGradient = 0;
        m_state.batchAccumulation = 0;
    }
}

void AdamOptimizer::step(Float gradient) {
    ++m_state.iter;

    Float actualLearningRate = m_hparams.learningRate * std::sqrt(1 - std::pow(m_hparams.beta2, m_state.iter)) / (1 - std::pow(m_hparams.beta1, m_state.iter));
    m_state.firstMoment = m_hparams.beta1 * m_state.firstMoment + (1 - m_hparams.beta1) * gradient;
    m_state.secondMoment = m_hparams.beta2 * m_state.secondMoment + (1 - m_hparams.beta2) * gradient * gradient;
    m_state.variable -= actualLearningRate * m_state.firstMoment / (std::sqrt(m_state.secondMoment) + m_hparams.epsilon);

    // Clamp the variable to the range [-20, 20] as a safeguard to avoid numerical instability:
    // since the sigmoid involves the exponential of the variable, value of -20 or 20 already yield
    // in *extremely* small and large results that are pretty much never necessary in practice.
    m_state.variable = std::min(std::max(m_state.variable, -20.0f), 20.0f);
}

Float AdamOptimizer::variable() const {
    return m_state.variable;
}

/* ============ */
/* QuadTreeNode */
/* ============ */

QuadTreeNode::QuadTreeNode() {
    m_children = {};
    for (size_t i = 0; i < m_sum.size(); ++i) {
        m_sum[i].store(0, std::memory_order_relaxed);
    }
}

QuadTreeNode::QuadTreeNode(const QuadTreeNode& arg) {
    copyFrom(arg);
}

QuadTreeNode& QuadTreeNode::operator=(const QuadTreeNode& arg) {
    copyFrom(arg);
    return *this;
}

void QuadTreeNode::setSum(int index, Float val) {
    m_sum[index].store(val, std::memory_order_relaxed);
}

void QuadTreeNode::setSum(Float val) {
    for (int i = 0; i < 4; ++i) {
        setSum(i, val);
    }
}

Float QuadTreeNode::sum(int index) const {
    return m_sum[index].load(std::memory_order_relaxed);
}

void QuadTreeNode::setChild(int idx, uint16_t val) {
    m_children[idx] = val;
}

uint16_t QuadTreeNode::child(int idx) const {
    return m_children[idx];
}

void QuadTreeNode::copyFrom(const QuadTreeNode& arg) {
    for (int i = 0; i < 4; ++i) {
        setSum(i, arg.sum(i));
        m_children[i] = arg.m_children[i];
    }
}

int QuadTreeNode::childIndex(Point2& p) const {
    int res = 0;
    for (int i = 0; i < Point2::dim; ++i) {
        if (p[i] < 0.5f) {
            p[i] *= 2;
        } else {
            p[i] = (p[i] - 0.5f) * 2;
            res |= 1 << i;
        }
    }

    return res;
}

Float QuadTreeNode::eval(Point2& p, const std::vector<QuadTreeNode>& nodes) const {
    SAssert(p.x >= 0 && p.x <= 1 && p.y >= 0 && p.y <= 1);
    
    const int index = childIndex(p);
    if (isLeaf(index)) {
        return 4 * sum(index);
    } else {
        return 4 * nodes[child(index)].eval(p, nodes);
    }
}

Float QuadTreeNode::pdf(Point2& p, const std::vector<QuadTreeNode>& nodes) const {
    SAssert(p.x >= 0 && p.x <= 1 && p.y >= 0 && p.y <= 1);
    const int index = childIndex(p);
    if (!(sum(index) > 0)) {
        return 0;
    }

    const Float factor = 4 * sum(index) / (sum(0) + sum(1) + sum(2) + sum(3));
    if (isLeaf(index)) {
        return factor;
    } else {
        return factor * nodes[child(index)].pdf(p, nodes);
    }
}

int QuadTreeNode::depthAt(Point2& p, const std::vector<QuadTreeNode>& nodes) const {
    SAssert(p.x >= 0 && p.x <= 1 && p.y >= 0 && p.y <= 1);

    const int index = childIndex(p);
    if (isLeaf(index)) {
        return 1;
    } else {
        return 1 + nodes[child(index)].depthAt(p, nodes);
    }
}

Point2 QuadTreeNode::sample(const std::vector<QuadTreeNode>& nodes) const {
    int index = 0;

    Float topLeft = sum(0);
    Float topRight = sum(1);
    Float partial = topLeft + sum(2);
    Float total = partial + topRight + sum(3);

    // Should only happen when there are numerical instabilities.
    if (!(total > 0.0f)) {
        return DirectionalTree::random->next2D();
    }

    Float boundary = partial / total;
    Point2 origin = Point2{0.0f, 0.0f};

    Float sample = DirectionalTree::random->next1D();

    if (sample < boundary) {
        SAssert(partial > 0);
        sample /= boundary;
        boundary = topLeft / partial;
    } else {
        partial = total - partial;
        SAssert(partial > 0);
        origin.x = 0.5f;
        sample = (sample - boundary) / (1.0f - boundary);
        boundary = topRight / partial;
        index |= 1 << 0;
    }

    if (sample < boundary) {
        sample /= boundary;
    } else {
        origin.y = 0.5f;
        sample = (sample - boundary) / (1.0f - boundary);
        index |= 1 << 1;
    }

    if (isLeaf(index)) {
        return origin + 0.5f * DirectionalTree::random->next2D();
    } else {
        return origin + 0.5f * nodes[child(index)].sample(nodes);
    }
}

void QuadTreeNode::record(Point2& p, Float irradiance, std::vector<QuadTreeNode>& nodes) {
    SAssert(p.x >= 0 && p.x <= 1 && p.y >= 0 && p.y <= 1);
    int index = childIndex(p);

    if (isLeaf(index)) {
        addToAtomicFloat(m_sum[index], irradiance);
    } else {
        nodes[child(index)].record(p, irradiance, nodes);
    }
}

void QuadTreeNode::record(const Point2& origin, Float size, Point2 nodeOrigin, Float nodeSize, Float value, std::vector<QuadTreeNode>& nodes) {
    Float childSize = nodeSize / 2;
    for (int i = 0; i < 4; ++i) {
        Point2 childOrigin = nodeOrigin;
        if (i & 1) { childOrigin[0] += childSize; }
        if (i & 2) { childOrigin[1] += childSize; }

        Float w = computeOverlappingArea(origin, origin + Point2(size), childOrigin, childOrigin + Point2(childSize));
        if (w > 0.0f) {
            if (isLeaf(i)) {
                addToAtomicFloat(m_sum[i], value * w);
            } else {
                nodes[child(i)].record(origin, size, childOrigin, childSize, value, nodes);
            }
        }
    }
}

Float QuadTreeNode::computeOverlappingArea(const Point2& min1, const Point2& max1, const Point2& min2, const Point2& max2) {
    Float lengths[2];
    for (int i = 0; i < 2; ++i) {
        lengths[i] = std::max(std::min(max1[i], max2[i]) - std::max(min1[i], min2[i]), 0.0f);
    }
    return lengths[0] * lengths[1];
}

bool QuadTreeNode::isLeaf(int index) const {
    return child(index) == 0;
}

void QuadTreeNode::build(std::vector<QuadTreeNode>& nodes) {
    for (int i = 0; i < 4; ++i) {
        // During sampling, all irradiance estimates are accumulated in
        // the leaves, so the leaves are built by definition.
        if (isLeaf(i)) continue;

        QuadTreeNode& c = nodes[child(i)];

        // Recursively build each child such that their sum becomes valid...
        c.build(nodes);

        // ...then sum up the children's sums.
        Float sum = 0;
        for (int j = 0; j < 4; ++j) {
            sum += c.sum(j);
        }
        setSum(i, sum);
    }
}

/* ============= */
/* InternalDTree */
/* ============= */

InternalDTree::InternalDTree() {
    m_atomic.sum.store(0, std::memory_order_relaxed);
    m_maxDepth = 0;
    m_nodes.emplace_back();
    m_nodes.front().setSum(0.0f);
}

const QuadTreeNode& InternalDTree::node(size_t i) const {
    return m_nodes[i];
}

Float InternalDTree::mean() const {
    if (m_atomic.statisticalWeight == 0) return 0;

    const Float factor = 1 / (M_PI * 4 * m_atomic.statisticalWeight);
    return factor * m_atomic.sum;
}

void InternalDTree::recordIrradiance(Point2 p, Float irradiance, Float statisticalWeight, EDirectionalFilter directionalFilter) {
    if (std::isfinite(statisticalWeight) && statisticalWeight > 0) {
        addToAtomicFloat(m_atomic.statisticalWeight, statisticalWeight);

        if (std::isfinite(irradiance) && irradiance > 0) {
            if (directionalFilter == EDirectionalFilter::ENearest) {
                m_nodes[0].record(p, irradiance * statisticalWeight, m_nodes);
            } else {
                int depth = depthAt(p);
                Float size = std::pow(0.5f, depth);

                Point2 origin = p;
                origin.x -= size / 2;
                origin.y -= size / 2;
                m_nodes[0].record(origin, size, Point2(0.0f), 1.0f, irradiance * statisticalWeight / (size * size), m_nodes);
            }
        }
    }
}

Float InternalDTree::pdf(Point2 p) const {
    if (!(mean() > 0)) {
        return 1 / (4 * M_PI);
    }

    return m_nodes[0].pdf(p, m_nodes) / (4 * M_PI);
}

int InternalDTree::depthAt(Point2 p) const {
    return m_nodes[0].depthAt(p, m_nodes);
}

int InternalDTree::depth() const {
    return m_maxDepth;
}

Point2 InternalDTree::sample() const {
    if (!(mean() > 0)) {
        return DirectionalTree::random->next2D();
    }

    Point2 res = m_nodes[0].sample(m_nodes);

    res.x = math::clamp(res.x, 0.0f, 1.0f);
    res.y = math::clamp(res.y, 0.0f, 1.0f);

    return res;
}

size_t InternalDTree::numNodes() const {
    return m_nodes.size();
}

Float InternalDTree::statisticalWeight() const {
    return m_atomic.statisticalWeight;
}

void InternalDTree::setStatisticalWeight(Float statisticalWeight) {
    m_atomic.statisticalWeight = statisticalWeight;
}

void InternalDTree::reset(const InternalDTree& previousDTree, int newMaxDepth, Float subdivisionThreshold) {
    m_atomic = Atomic{};
    m_maxDepth = 0;
    m_nodes.clear();
    m_nodes.emplace_back();

    struct StackNode {
        size_t nodeIndex;
        size_t otherNodeIndex;
        const InternalDTree* otherDTree;
        int depth;
    };

    std::stack<StackNode> nodeIndices;
    nodeIndices.push({0, 0, &previousDTree, 1});

    const Float total = previousDTree.m_atomic.sum;
    
    // Create the topology of the new DTree to be the refined version
    // of the previous DTree. Subdivision is recursive if enough energy is there.
    while (!nodeIndices.empty()) {
        StackNode sNode = nodeIndices.top();
        nodeIndices.pop();

        m_maxDepth = std::max(m_maxDepth, sNode.depth);

        for (int i = 0; i < 4; ++i) {
            const QuadTreeNode& otherNode = sNode.otherDTree->m_nodes[sNode.otherNodeIndex];
            const Float fraction = total > 0 ? (otherNode.sum(i) / total) : std::pow(0.25f, sNode.depth);
            SAssert(fraction <= 1.0f + Epsilon);

            if (sNode.depth < newMaxDepth && fraction > subdivisionThreshold) {
                if (!otherNode.isLeaf(i)) {
                    SAssert(sNode.otherDTree == &previousDTree);
                    nodeIndices.push({m_nodes.size(), otherNode.child(i), &previousDTree, sNode.depth + 1});
                } else {
                    nodeIndices.push({m_nodes.size(), m_nodes.size(), this, sNode.depth + 1});
                }

                m_nodes[sNode.nodeIndex].setChild(i, static_cast<uint16_t>(m_nodes.size()));
                m_nodes.emplace_back();
                m_nodes.back().setSum(otherNode.sum(i) / 4);

                if (m_nodes.size() > std::numeric_limits<uint16_t>::max()) {
                    SLog(EWarn, "Maximum children count hit.");
                    nodeIndices = std::stack<StackNode>();
                    break;
                }
            }
        }
    }

    // Uncomment once memory becomes an issue.
    //m_nodes.shrink_to_fit();

    for (auto& node : m_nodes) {
        node.setSum(0);
    }
}

size_t InternalDTree::approxMemoryFootprint() const {
    return m_nodes.capacity() * sizeof(QuadTreeNode) + sizeof(*this);
}

void InternalDTree::build() {
    auto& root = m_nodes[0];

    // Build the quadtree recursively, starting from its root.
    root.build(m_nodes);

    // Ensure that the overall sum of irradiance estimates equals
    // the sum of irradiance estimates found in the quadtree.
    Float sum = 0;
    for (int i = 0; i < 4; ++i) {
        sum += root.sum(i);
    }
    m_atomic.sum.store(sum);
}

/* =============== */
/* DirectionalTree */
/* =============== */

ref<RandomGen> DirectionalTree::random = new RandomGen();

void DirectionalTree::record(const DTreeRecord& rec, EDirectionalFilter directionalFilter, EBsdfSamplingFractionLoss bsdfSamplingFractionLoss) {
    if (!rec.isDelta) {
        Float irradiance = rec.radiance / rec.woPdf;

        Vector dir = rec.d;
        Point2 p{0, 0};

        if (std::isfinite(dir.x) && std::isfinite(dir.y) && std::isfinite(dir.z)) {
            const Float cosTheta = std::min(std::max(dir.z, -1.0f), 1.0f);
            Float phi = std::atan2(dir.y, dir.x);
            while (phi < 0)
                phi += 2.0 * M_PI;

            p = Point2{(cosTheta + 1) / 2, phi / (2 * M_PI)};
        }

        building.recordIrradiance(p, irradiance, rec.statisticalWeight, directionalFilter);
    }

    if (bsdfSamplingFractionLoss != EBsdfSamplingFractionLoss::ENone && rec.product > 0) {
        optimizeBsdfSamplingFraction(rec, bsdfSamplingFractionLoss == EBsdfSamplingFractionLoss::EKL ? 1.0f : 2.0f);
    }
}

void DirectionalTree::build() {
    building.build();
    sampling = building;
}

void DirectionalTree::reset(int maxDepth, Float subdivisionThreshold) {
    building.reset(sampling, maxDepth, subdivisionThreshold);
}

Vector DirectionalTree::sample() const {
    Point2 sample = sampling.sample();

    const Float cosTheta = 2 * sample.x - 1;
    const Float phi = 2 * M_PI * sample.y;

    const Float sinTheta = sqrt(1 - cosTheta * cosTheta);
    Float sinPhi, cosPhi;
    math::sincos(phi, &sinPhi, &cosPhi);

    return {sinTheta * cosPhi, sinTheta * sinPhi, cosTheta};
}

Float DirectionalTree::pdf(const Vector& dir) const {
    Point2 p{0, 0};

    if (std::isfinite(dir.x) && std::isfinite(dir.y) && std::isfinite(dir.z)) {
        const Float cosTheta = std::min(std::max(dir.z, -1.0f), 1.0f);
        Float phi = std::atan2(dir.y, dir.x);
        while (phi < 0)
            phi += 2.0 * M_PI;

        p = Point2{(cosTheta + 1) / 2, phi / (2 * M_PI)};
    }

    return sampling.pdf(p);
}

Float DirectionalTree::pdf(const Point2& uv) const {
    return sampling.pdf(uv);
}

int DirectionalTree::depth() const {
    return sampling.depth();
}

size_t DirectionalTree::numNodes() const {
    return sampling.numNodes();
}

Float DirectionalTree::meanRadiance() const {
    return sampling.mean();
}

Float DirectionalTree::statisticalWeight() const {
    return sampling.statisticalWeight();
}

Float DirectionalTree::statisticalWeightBuilding() const {
    return building.statisticalWeight();
}

void DirectionalTree::setStatisticalWeightBuilding(Float statisticalWeight) {
    building.setStatisticalWeight(statisticalWeight);
}

size_t DirectionalTree::approxMemoryFootprint() const {
    return building.approxMemoryFootprint() + sampling.approxMemoryFootprint();
}

void DirectionalTree::optimizeBsdfSamplingFraction(const DTreeRecord& rec, Float ratioPower) {
    m_lock.lock();

    // GRADIENT COMPUTATION
    Float variable = bsdfSamplingFractionOptimizer.variable();
    Float samplingFraction = bsdfSamplingFraction(variable);

    // Loss gradient w.r.t. sampling fraction
    Float mixPdf = samplingFraction * rec.bsdfPdf + (1 - samplingFraction) * rec.dTreePdf;
    Float ratio = std::pow(rec.product / mixPdf, ratioPower);
    Float dLoss_dSamplingFraction = -ratio / rec.woPdf * (rec.bsdfPdf - rec.dTreePdf);

    // Chain rule to get loss gradient w.r.t. trainable variable
    Float dLoss_dVariable = dLoss_dSamplingFraction * dBsdfSamplingFraction_dVariable(variable);

    // We want some regularization such that our parameter does not become too big.
    // We use l2 regularization, resulting in the following linear gradient.
    Float l2RegGradient = 0.01f * variable;

    Float lossGradient = l2RegGradient + dLoss_dVariable;

    // ADAM GRADIENT DESCENT
    bsdfSamplingFractionOptimizer.append(lossGradient, rec.statisticalWeight);

    m_lock.unlock();
}

void DirectionalTree::dump(BlobWriter& blob, const Point& p, const Vector& size) const {
    blob
        << (float)p.x << (float)p.y << (float)p.z
        << (float)size.x << (float)size.y << (float)size.z
        << (float)sampling.mean() << (uint64_t)sampling.statisticalWeight() << (uint64_t)sampling.numNodes();

    for (size_t i = 0; i < sampling.numNodes(); ++i) {
        const auto& node = sampling.node(i);
        for (int j = 0; j < 4; ++j) {
            blob << (float)node.sum(j) << (uint16_t)node.child(j);
        }
    }
}

void DirectionalTree::construct(DSArguments& init_data)
{
    return;
}

void DirectionalTree::preprocess()
{
    return;
}

void DirectionalTree::store(std::vector<Sample>& samples)
{
    // As Müller et al.'s D-Trees require multiple iterations, we first store our sample contingent
    // into a vector based on a geometric series, then use these samples in the postprocessing step
    // to build the actual tree.

    auto total_samples = samples.size();
    for (size_t s_i = 0; s_i < total_samples; ++s_i)
    {
        Sample sample = samples.at(s_i);

        Point2 p(
            0.5 * (std::cos(sample.theta) + 1),
            INV_TWOPI * sample.phi
        );

        sample.theta = p.x;
        sample.phi = p.y;

        this->l_sample_storage.push_back(sample);
    }
}

void DirectionalTree::postprocess()
{
    // Müller et al. states that the geometric series uses twice as many samples as in the previous
    // iteration. We therefore iterate over all stored samples, building and resetting when we hit a threshold.
    // We start with 1 sample and go from there.

    size_t t = 1;

    auto total_samples = this->l_sample_storage.size();
    for (size_t s_i = 0; s_i < total_samples; ++s_i)
    {
        Sample sample = this->l_sample_storage.at(s_i);
        building.recordIrradiance(Point2(sample.theta, sample.phi), sample.value, INV_FOURPI, EDirectionalFilter::ENearest);

        //optimizeBsdfSamplingFraction()

        if ((s_i == t - 1) && (total_samples - t >= t))
        {
            build();
            reset(20, 0.01f);
            t *= 2;
        }
    }

    build();

    //std::cout << "depth: " << sampling.depth() << std::endl;
    //std::cout << "nodes: " << sampling.numNodes() << std::endl;
    //std::cout << "mean: " << sampling.mean() << std::endl;
}

Sample DirectionalTree::sample(Point2& pos)
{
    // We ignore the sample we pass in and instead use the internal sample/) function.
    Point2 coords = sampling.sample();
    
    Float theta = std::acos(2 * coords.x - 1);
    Float phi = 2 * M_PI * coords.y;

    Float pdf = sampling.pdf(
        Point2(
            0.5 * (std::cos(theta) + 1),
            INV_TWOPI * phi
        )
    );

    Sample sample = {
        .value = pdf,
        .phi = phi,
        .theta = theta
    };
    return sample;
}

void DirectionalTree::wipe()
{
    building = InternalDTree();
    this->l_sample_storage.clear();
}

DSType DirectionalTree::type()
{
    return DSType::DS_DTree;
}

MTS_NAMESPACE_END