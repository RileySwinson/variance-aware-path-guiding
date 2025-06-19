#pragma once

#if !defined(__DSCOMPARE_UTIL_STATTRAK_H_)
#define __DSCOMPARE_UTIL_STATTRAK_H_

#include <ds-compare/util/_definitions.h>

MTS_NAMESPACE_BEGIN

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
        const auto size = this->m_data.size();
        for (int i = 0; i < size; ++i)
        {
            this->m_data[static_cast<DSType>(i)] = StatData();
        }
    }

    /// Completely clears the tracker.
    void clear()
    {
        this->m_data.clear();
    }

    void timer_start(const std::string& identifier = "")
    {
        /* Check if there's even an entry for the currently active ds, if not, create one */
        auto it = this->m_data.find(this->m_active);
        if (it == this->m_data.end())
        {
            this->m_data[this->m_active];
        }
        
        StatData& data = this->m_data.at(this->m_active);
        data.m_times.emplace(identifier, std::vector<TimeMeasure>());

        TimeMeasure& measure = get_first_available(data.m_times.at(identifier), true);
        measure.start = std::chrono::steady_clock::now();
    }

    void timer_end(const std::string& identifier = "")
    {
        auto end_measure = std::chrono::steady_clock::now();

        StatData& data = this->m_data.at(this->m_active);

        auto it = data.m_times.find(identifier);
        if (it == data.m_times.end())
        {
            SLog(EError, "Identifier '%s' could not be found.", identifier.c_str());
        }

        TimeMeasure& measure = get_first_available(it->second, false);
        measure.end = end_measure;
        measure.finished = true;
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

        /* Initialize header with errors + time measurements */
        std::string header = ",MD,RMSE,MSE,MAE,Memory,";
        const auto& f_times = this->m_data.begin()->second.m_times;

        std::vector<std::string> time_headers;
        for (const auto& value : f_times)
        {
            const std::string time_header = value.first;
            header += time_header + " (s),";
            time_headers.push_back(time_header);
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

    TimeMeasure& get_first_available(std::vector<TimeMeasure>& slots, bool rec_start)
    {
        for (auto& slot : slots)
        {
            if (slot.finished) continue;
            if (rec_start && slot.start) SLog(EError, "Attempting to start a timer before finishing the previous segment.");

            return slot;
        }

        if (!rec_start) SLog(EError, "Attempting to stop a timer although there is no active timer running.");

        slots.push_back(TimeMeasure());
        return slots.back();
    }

    Float get_duration_helper(const TimesMap& times, const std::string& identifier, const int slot)
    {
        auto it = times.find(identifier);
        if (it == times.end())
            SLog(EError, "Identifier '%s' could not be found.", identifier.c_str());

        auto tm = it->second;

        if (slot < -1 || slot >= (int) tm.size())
            SLog(EError, "Timing slot index is out of bounds in StatTrak::get_duration");

        Float total_duration = 0;
        for (int i = 0; i < tm.size(); ++i)
        {
            if (slot != -1 && slot != i) continue;

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