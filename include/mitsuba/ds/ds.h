#pragma once

#if !defined(__MITSUBA_DS_DS_H_)
#define __MITSUBA_DS_DS_H_

#include <mitsuba/mitsuba.h>
#include <memory>
#include <map>
#include <functional>

MTS_NAMESPACE_BEGIN

/**
 * \brief Enum values for the respective data structures used in the dscompare plugin.
 * 
 * This enum provides identifiers required for both storage and identification in the \ref DSCluster.
 * If a value is missing, the missing value may be added by the user.
 */
enum MTS_EXPORT_CORE DSType {
    DS_GaussianMixture,
    DS_TileCoding,
    DS_SphericalHarmonics,
    DS_DTree,
    DS_Invalid
};

/**
 * Storage struct for the various parameters a data structure must need. This is passed
 * into each construct(), where data structures can then individually fetch the needed information for initialization.
 * If a value is missing, the missing value may be added by the user.
 */
struct MTS_EXPORT_CORE DSInitData {
	int sh_bands = 3;
};

struct MTS_EXPORT_CORE Sample {
    float luminance = 0.0f;
    
    float phi = 0.0f;
    float theta = 0.0f;
};

/**
 * \brief Abstract base class for all data structures.
 * 
 * As the approach of the comparison framework itself is based on a plug-and-play plugin system, each new
 * data structure must extend this class and implement the pure virtual methods to be usable within the comparison framework.
 */
struct MTS_EXPORT_CORE DataStructure {
    virtual ~DataStructure() { }

    /// Calls all relevant functions and constructs the data structure in such a way that it is ready-to-use for data storage.
    virtual DataStructure* construct(DSInitData& init_data) = 0;

    /// Stores a sample into the data structure.
    virtual void store(Sample& sample) = 0;

    /// Obtain a sample from the underlying approximation that is stored in the data structure.
    virtual Sample sample(Point2& pos) = 0;

    /// Clears the entire data structure such that it is back to its initial, empty state.
    virtual void wipe() = 0;

    /// \brief Returns the type of the data structure, see \ref DSType.
    virtual DSType type() = 0;
};

/**
 * \brief Manager class for data structure storage and handling.
 * 
 * A singleton responsible 
 */
class MTS_EXPORT_CORE DSCluster {
public:
    typedef std::map<DSType, DataStructure*>::iterator DSIter;

    /// Returns an instance of the data structure cluster.
    static inline DSCluster& get()
    {
        static DSCluster instance;
        return instance;
    };

    /// Attaches a data structure instance to the cluster.
    bool attach(DataStructure* ds)
    {
        return m_map.emplace(ds->type(), ds).second;
    }

    /**
     * \brief Attaches a data structure instance to the cluster.
     * 
     * Please note that this is an overload that should only be used if the user intends to overwrite
     * or manually specify the DSType with which the data structure is stored and can be referenced with.
     * Otherwise, please refer to \ref attach(DataStructure*)
     */
    bool attach(const DSType type, DataStructure* ds)
    {
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

#endif /* __MITSUBA_DS_DS_H_ */