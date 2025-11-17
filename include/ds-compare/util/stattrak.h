#pragma once

#if !defined(__DSCOMPARE_UTIL_STATTRAK_H_)
#define __DSCOMPARE_UTIL_STATTRAK_H_

#include <ds-compare/util/_definitions.h>

MTS_NAMESPACE_BEGIN

#if !defined(STATTRAK_FUNCTION_TIMER)
    #define STATTRAK_FUNCTION_TIMER(name) \
       auto stattrak_timer_##__LINE__ = StatTrak::get().scoped_timer(name)
#endif

#if !defined(STATTRAK_BLOCK_TIMER)
    #define STATTRAK_BLOCK_TIMER(name) \
        if (auto stattrak_timer_##__LINE__ = StatTrak::get().scoped_timer(name))
#endif

enum DS_COMPARE DSType : int;

enum DS_COMPARE MeasureMetric : int {
    MD,
    RMSE,
    MSE,
    MAE,
    Memory
};

struct TimeMeasure {
    boost::optional<std::chrono::steady_clock::time_point> start;
    boost::optional<std::chrono::steady_clock::time_point> end;
    bool finished = false;
};

struct DS_COMPARE StatTrak {
    static inline StatTrak& get()
    {
        static StatTrak instance;
        return instance;
    }

    void follow(DSType ds)
    {
        this->m_active = ds;
    }

    void store(MeasureMetric metric, Float value)
    {
        auto it = this->m_data.find(this->m_active);
        if (it == this->m_data.end())
        {
            this->m_data[this->m_active];
        }
        
        this->m_data.at(this->m_active).m_errors.emplace(metric, value);
    }

    /// Resets all entries for each data structure within the tracker back to 0.
    void reset()
    {
        for (auto& entry : this->m_data)
        {
            entry.second = StatData();
        }
    }

    /// Completely clears the tracker.
    void clear()
    {
        this->m_data.clear();
    }

    /// Specify which identifiers should be used in the final output. If not specified, all will be used.
    template<typename... Args>
    void set_valid_timers(Args&&... args)
    {
        (void)std::initializer_list<int>{
            (this->whitelist.push_back(std::forward<Args>(args)), 0)...
        };
    }

    /// Obtain the duration between two measured points specified by the identifier for the currently assigned data structure.
    Float get_duration(const std::string& identifier = "", const int slot = -1)
    {
        auto it = this->m_data.find(this->m_active);
        if (it == this->m_data.end())
        {
            SLog(EError, "Data structure could not be found.");
        }

        auto& times = this->m_data.at(this->m_active).m_times;
        return get_duration_helper(times, identifier, slot);
    }

    /// Obtain the duration between two measured points specified by the identifier for the specified data structure.
    Float get_duration(const DSType ds, const std::string& identifier = "", const int slot = -1)
    {
        auto data_it = this->m_data.find(ds);
        if (data_it == this->m_data.end())
        {
            SLog(EError, "Data structure could not be found.");
        }

        auto& times = data_it->second.m_times;
        return get_duration_helper(times, identifier, slot);
    }

    /// Obtain the error for the currently assigned data structure.
    Float get_error(const MeasureMetric metric)
    {
        auto it = this->m_data.find(this->m_active);
        if (it == this->m_data.end())
        {
            SLog(EError, "Data structure could not be found.");
        }

        auto& errors = this->m_data.at(this->m_active).m_errors;
        return get_error_helper(errors, metric);
    }

    /// Obtain the error for the specified data structure.
    Float get_error(const DSType ds, const MeasureMetric metric)
    {
        auto data_it = this->m_data.find(ds);
        if (data_it == this->m_data.end())
        {
            SLog(EError, "Data structure could not be found.");
        }

        auto& errors = data_it->second.m_errors;
        return get_error_helper(errors, metric);
    }

    /// Write the content of the tracker to an output file.
    void write(const std::string& path)
    {
        if (this->m_data.empty()) return;

        std::ofstream output;
        output.open(path, std::ios::out);

        /* Accumulate timer headers */
        std::set<std::string> time_headers;
        for (const auto& entry : this->m_data)
        {
            for (const auto& time_pair : entry.second.m_times)
            {
                time_headers.insert(time_pair.first);
            }
        }

        /* Initialize header with errors + time measurements */
        std::string header = ",MD,RMSE,MSE,MAE,Memory,";
        for (const auto& time_header : time_headers)
        {
            header += time_header + " (s),";
        }
        header += "\n";

        output << header;

        /* Iterate over all entries and accumulate the stored data into an output string */
        std::string res = "";
        for (const auto& entry : this->m_data)
        {
            DSType type = entry.first;
            res += std::to_string(type) + ",";

            const auto& errors = entry.second.m_errors;
            const auto& times = entry.second.m_times;

            auto add_measure = [&](MeasureMetric m) { res += (errors.find(m) != errors.end()) ? std::to_string(get_error(type, m)) + "," : ","; };

            // Errors
            add_measure(MD);
            add_measure(RMSE);
            add_measure(MSE);
            add_measure(MAE);

            // Memory
            add_measure(Memory);

            // Time
            for (const auto& header : time_headers)
            {
                res += (times.find(header) != times.end()) ? std::to_string(get_duration(type, header)) + "," : ",";
            }
            res += "\n";
        }

        output << res;
        output.close();
    }


    class ScopedTimer {
    public:
        ScopedTimer(StatTrak& parent, const std::string& identifier) : m_parent(parent), m_identifier(identifier)
        {
            this->m_index = this->m_parent.timer_start(this->m_identifier);
        }

        ~ScopedTimer()
        {
            this->m_parent.timer_end(this->m_identifier, this->m_index);
        }

        ScopedTimer(const ScopedTimer&) = delete;
        ScopedTimer& operator=(const ScopedTimer&) = delete;

        ScopedTimer(ScopedTimer&&) = default;
        ScopedTimer& operator=(ScopedTimer&&) = default;

        explicit operator bool() const { return true; }

    private:
        StatTrak& m_parent;
        const std::string m_identifier;
        size_t m_index;
    };

    ScopedTimer scoped_timer(const std::string& identifier)
    {
        return ScopedTimer(*this, identifier);
    }

    StatTrak(StatTrak const&)       = delete;
    void operator=(StatTrak const&) = delete;
private:
    using ErrorMap = std::map<MeasureMetric, Float>;
    using TimesMap = std::map<std::string, std::vector<TimeMeasure>>;

    struct StatData {
        ErrorMap m_errors;
        TimesMap m_times;
    };

    DSType m_active;
    std::map<DSType, StatData> m_data;
    std::vector<std::string> whitelist;

    size_t timer_start(const std::string& identifier = "")
    {
        /* Cancel if the identifier isn't whitelisted */
        if (std::find(this->whitelist.begin(), this->whitelist.end(), identifier) == this->whitelist.end())
        {
            return 0;
        }

        /* Check if there's even an entry for the currently active ds, if not, create one */
        auto it = this->m_data.find(this->m_active);
        if (it == this->m_data.end())
        {
            this->m_data[this->m_active];
        }

        std::vector<TimeMeasure>& slots = this->m_data.at(this->m_active).m_times[identifier];
        slots.push_back(TimeMeasure());

        slots.back().start = std::chrono::steady_clock::now();

        return slots.size() - 1;
    }

    void timer_end(const std::string& identifier = "", size_t slot_index = 0)
    {
        /* Cancel if the identifier isn't whitelisted */
        if (std::find(this->whitelist.begin(), this->whitelist.end(), identifier) == this->whitelist.end())
        {
            return;
        }

        auto end_measure = std::chrono::steady_clock::now();

        StatData& data = this->m_data.at(this->m_active);

        auto it = data.m_times.find(identifier);
        if (it == data.m_times.end())
        {
            SLog(EError, "Identifier '%s' could not be found.", identifier.c_str());
            return;
        }

        std::vector<TimeMeasure>& slots = it->second;
        if (slot_index >= slots.size())
        {
            SLog(EError, "Timer slot index %zu is out of bounds for identifier '%s'.", slot_index, identifier.c_str());
            return;
        }

        TimeMeasure& measure = slots[slot_index];
        if (measure.finished)
        {
            SLog(EError, "Timer for '%s' at slot %zu was already finished.", identifier.c_str(), slot_index);
            return;
        }

        measure.end = end_measure;
        measure.finished = true;
    }

    Float get_duration_helper(const TimesMap& times, const std::string& identifier, const int slot)
    {
        auto it = times.find(identifier);
        if (it == times.end())
            SLog(EError, "Identifier '%s' could not be found.", identifier.c_str());

        auto& tm = it->second;

        if (slot < -1 || slot >= (int) tm.size())
            SLog(EError, "Timing slot index is out of bounds in StatTrak::get_duration");

        Float total_duration = 0;
        for (size_t i = 0; i < tm.size(); ++i)
        {
            if (slot != -1 && (size_t) slot != i) continue;

            auto& measure = tm.at(i);
            if (!measure.finished) continue;

            auto diff = measure.end.get() - measure.start.get();
            total_duration += std::chrono::duration_cast<std::chrono::duration<Float>>(diff).count();
        }

        return total_duration;
    }

    Float get_error_helper(const ErrorMap& errors, const MeasureMetric metric)
    {
        auto it = errors.find(metric);
        if (it == errors.end())
            SLog(EError, "Error metric could not be found.");

        return it->second;
    }

    StatTrak() = default;
    ~StatTrak()
    {
        clear();
    }
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_UTIL_STATTRAK_H_ */