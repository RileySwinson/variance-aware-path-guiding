#pragma once

#if !defined(__DSCOMPARE_DS_H_)
#define __DSCOMPARE_DS_H_

#include <ds-compare/util/_util.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Enum values for the respective data structures used in the dscompare plugin.
 *
 * This enum provides identifiers required for both storage and identification in the DSCluster.
 * If a value is missing, the missing value may be added by the user.
 */
enum DS_COMPARE DSType : int {
    DS_Invalid = -1,
    DS_Unidirectional,
    DS_SphericalHarmonics,
    DS_DTree,
    DS_VMFMixture,
    DS_TileCoding,
    DS_BinaryTileCoding
};

/**
 * Storage struct for the various parameters a data structure must need. This is passed
 * into each construct(), where data structures can then individually fetch the needed information for initialization.
 * If a value is missing, the missing value may be added by the user.
 */
struct DS_COMPARE DSArguments {
    struct Comparer {
        std::string path = "./data/tests/envmaps/";
        std::string result_path = "./data/results";
        Sample::Mode mode = Sample::Mode::Sphere;
        uint32_t samples_learning = 16384; // 2^14
        uint32_t samples_guiding = 65536; // 2^16
        std::vector<int> blacklist;
        bool normalize = false;
        bool visualize = true;
        EnvironmentMap::VisualizationMode vis_mode = EnvironmentMap::VisualizationMode::Mono;
    };

    struct Noise {
        bool envmap = false;
        bool samples = false;
    };

    struct SH {
        int bands = 7;
        int depth = 12;
        bool use_offset = true;
    };
    
    struct DTree {
        DTreeParams::EBsdfSamplingFractionLoss frac_loss = DTreeParams::EBsdfSamplingFractionLoss::ENone;
        DTreeParams::EDirectionalFilter dir_filter = DTreeParams::EDirectionalFilter::ENearest;
        Float threshold = 0.01;
        int iterations = -1;
        int max_depth = 20;
    };

    struct TC {
        int tilings = 4;
        int tiles_x = 16;
        int tiles_y = 16;
    };

    struct BTC {
        int tilings = 4;
        int tiles_x = 1;
        int tiles_y = 1;
        int max_depth = 10;
        int eagerness = 0;
        TCParams::Transformation transformation_mode = TCParams::Transformation::Spherical;
    };

    struct VMF {
        uint32_t components = 16;
        bool use_ruppert = true;
    };

    Comparer comparer;
    Noise noise;
    SH sh;
    DTree dt;
    TC tc;
    BTC btc;
    VMF vmf;
};

/**
 * @brief Abstract base class for all data structures.
 * 
 * As the approach of the comparison framework itself is based on a plug-and-play plugin system, each new
 * data structure must extend this class and implement the pure virtual methods to be usable within the comparison framework.
 */
struct DS_COMPARE DataStructure {
    virtual ~DataStructure() { }

    /// Calls all relevant functions and initializes the data structure in such a way that it is ready-to-use for data storage.
    virtual void construct(DSArguments& init_data) = 0;

    /// Performs operations after construction but before storage, if necessary.
    virtual void preprocess() = 0;

    /// Stores a number of samples into the data structure.
    virtual void store(std::vector<Sample>& samples) = 0;

    /// Performs operations after storage but before sampling, if necessary.
    virtual void postprocess() = 0;

    /// Obtains a sample from the underlying approximation that is stored in the data structure.
    virtual Sample sample(Point2& pos) = 0;

    /// Evaluates the underlying sampling pdf at a given position in the domain [0, 1)^2.
    virtual Float eval(Point2& pos) = 0;

    /// Clears the entire data structure such that it is back to its initial, empty state.
    virtual void wipe() = 0;

    /// Returns the type of the data structure, see DSType.
    virtual DSType type() = 0;

    /// Returns the name of the data structure.
    virtual std::string name() = 0;

    /// Obtains the approximate memory footprint for this data structure (in bytes).
    virtual int memory() = 0;

    /* ==== Miscellaneous ==== */

    template<class Iterator>
    inline bool is_in(Iterator start, Iterator end)
    {
        if (start == end) return false;
        return (std::find(start, end, type()) != end);
    }
};

/**
 * @brief Manager class for data structure storage and handling.
 * 
 * A singleton responsible for performing operations on all data structures stored within it.
 */
class DS_COMPARE DSCluster {
public:
    typedef std::map<DSType, DataStructure*>::iterator DSIter;

    /// Returns an instance of the data structure cluster.
    static inline DSCluster& get()
    {
        static DSCluster instance;
        return instance;
    }

    /// Attaches a data structure instance to the cluster.
    bool attach(DataStructure* ds)
    {
        if (ds->type() == DSType::DS_Invalid) 
        {
            SLog(EWarn, "Attempting to attach data structure with type 'DS_Invalid' -- ignoring.");
            return false;
        }
        return m_map.emplace(ds->type(), ds).second;
    }

    /**
     * @brief Attaches a data structure instance to the cluster.
     * 
     * Please note that this is an overload that should only be used if the user intends to overwrite
     * or manually specify the DSType with which the data structure is stored and can be referenced with.
     * Otherwise, please refer to attach(DataStructure*)
     */
    bool attach(const DSType type, DataStructure* ds)
    {
        if (type == DSType::DS_Invalid || ds->type() == DSType::DS_Invalid) 
        {
            SLog(EWarn, "Attempting to attach data structure with type 'DS_Invalid' -- ignoring.");
            return false;
        }
        return m_map.emplace(type, ds).second;
    }

    /// Fetch the pointer to the data structure from the cluster given a DSType. If no data structure is found, a nullptr is returned.
    DataStructure* obtain(const DSType type) const
    {
        auto it = m_map.find(type);
        if (it == m_map.end()) 
        {
            SLog(ELogLevel::EWarn, ("Data structure (" + std::to_string(type) + ") unable to fetch. Make sure you have registered it correctly -- skipping.").c_str());
            return nullptr;
        }

        return it->second;
    }

    /// Applies a given input function F to all data structures within the cluster.
    void for_each(std::function<void(DataStructure*)> F)
    {
        for (DSIter it = this->m_map.begin(); it != this->m_map.end(); ++it)
        {
            F(it->second);
        } 
    }

    /// Returns the number of data structures stored within the cluster.
    int size()
    {
        return this->m_map.size();
    }

    /// Reverts the cluster back to its blank state.
    void clear()
    {
        for (DSIter it = this->m_map.begin(); it != this->m_map.end(); ++it)
        {
            delete it->second;
        }

        this->m_map.clear();
    }

    DSCluster(DSCluster const&)         = delete;
    void operator=(DSCluster const&)    = delete;
private:
    std::map<DSType, DataStructure*> m_map;

    DSCluster() = default;
    ~DSCluster()
    {
        this->clear();
    }
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_DS_H_ */