#pragma once

#if !defined(__DSCOMPARE_UTIL_STATTRAK_H_)
#define __DSCOMPARE_UTIL_STATTRAK_H_

#include <ds-compare/util/_definitions.h>

MTS_NAMESPACE_BEGIN

enum DS_COMPARE DSType : int;

enum DS_COMPARE ErrorMetric : int {
    RMSE,
    MSE,
    MAE
};

struct TimeMeasure {
    std::chrono::steady_clock::time_point start;
    std::chrono::steady_clock::time_point end;
};

struct DS_COMPARE StatTrak {
    static inline StatTrak& get()
    {
        static StatTrak instance;
        return instance;
    }

    void reserve(int capacity)
    {
        for (int i = 0; i < capacity; ++i)
        {
            this->m_data[static_cast<DSType>(i)];
        }
    }

    void follow(DSType ds)
    {
        this->m_active = ds;
    }

    void store(ErrorMetric ds, Float value)
    {
        auto it = this->m_data.find(this->m_active);
        if (it == this->m_data.end())
        {
            this->m_data[this->m_active];
        }
        
        this->m_data.at(this->m_active).m_errors.emplace(ds, value);
    }

    /// Resets all entries for each envmap within the tracker back to 0.
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

    void timer_start(const std::string& identifier)
    {
        auto it = this->m_data.find(this->m_active);
        if (it == this->m_data.end())
        {
            this->m_data[this->m_active];
        }

        StatData& data = this->m_data.at(this->m_active);
        data.m_times.emplace(identifier, TimeMeasure());
        data.m_times[identifier].start = std::chrono::steady_clock::now();
    }

    void timer_end(const std::string& identifier)
    {
        auto end_measure = std::chrono::steady_clock::now();

        StatData& data = this->m_data.at(this->m_active);

        auto it = data.m_times.find(identifier);
        if (it == data.m_times.end())
        {
            SLog(EError, "Identifier '%s' could not be found.", identifier.c_str());
        }

        it->second.end = end_measure;
    }

    /// Obtain the duration between two measured points specified by the identifier for the currently assigned data structure.
    Float get_duration(const std::string& identifier)
    {
        auto it = this->m_data.find(this->m_active);
        if (it == this->m_data.end())
        {
            SLog(EError, "Data structure could not be found.");
        }

        auto& times = this->m_data.at(this->m_active).m_times;
        return get_duration_helper(times, identifier).count();
    }

    /// Obtain the duration between two measured points specified by the identifier for the specified data structure.
    Float get_duration(const std::string& identifier, const DSType ds)
    {
        auto data_it = this->m_data.find(ds);
        if (data_it == this->m_data.end())
        {
            SLog(EError, "Data structure could not be found.");
        }

        auto& times = data_it->second.m_times;
        return get_duration_helper(times, identifier).count();
    }

    /// Obtain the error for the currently assigned data structure.
    Float get_error(const ErrorMetric metric)
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
    Float get_error(const ErrorMetric metric, const DSType ds)
    {
        auto data_it = this->m_data.find(ds);
        if (data_it == this->m_data.end())
        {
            SLog(EError, "Data structure could not be found.");
        }

        auto& errors = data_it->second.m_errors;
        return get_error_helper(errors, metric);
    }

    void write(const std::string& path)
    {
        if (this->m_data.empty()) return;

        std::ofstream output;
        output.open(path, std::ios::out);

        /* Initialize header with errors + time measurements */
        std::string header = ",RMSE,MSE,MAE,";
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

            // Errors
            res += (errors.find(RMSE) != errors.end()) ? std::to_string(get_error(RMSE, type)) + "," : ",";
            res += (errors.find(MSE) != errors.end())  ? std::to_string(get_error(MSE, type)) + ","  : ",";
            res += (errors.find(MAE) != errors.end())  ? std::to_string(get_error(MAE, type)) + ","  : ",";

            // Time
            for (const auto& header : time_headers)
            {
                res += (times.find(header) != times.end()) ? std::to_string(get_duration(header, type)) + "," : ",";
            }
            res += "\n";
        }

        output << res;
        output.close();
    }

    StatTrak(StatTrak const&)       = delete;
    void operator=(StatTrak const&) = delete;
private:
    using ErrorMap = std::map<ErrorMetric, Float>;
    using TimesMap = std::map<std::string, TimeMeasure>;

    struct StatData {
        ErrorMap m_errors;
        TimesMap m_times;
    };

    DSType m_active;
    std::map<DSType, StatData> m_data;

    std::chrono::duration<Float> get_duration_helper(const TimesMap& times, const std::string& identifier)
    {
        auto it = times.find(identifier);
        if (it == times.end())
            SLog(EError, "Identifier '%s' could not be found.", identifier.c_str());

        TimeMeasure tm = it->second;
        return std::chrono::duration_cast<std::chrono::duration<Float>>(tm.end - tm.start);
    }

    Float get_error_helper(const ErrorMap& errors, const ErrorMetric metric)
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