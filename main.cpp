/*
====================================================================================
 Stock Market Technical Analysis Tool (RSI, SVI, MACD)
====================================================================================
This C++ program performs technical analysis on stock market data using key indicators:
- **Relative Strength Index (RSI)**
- **Stochastic Volatility Index (SVI)**
- **Moving Average Convergence Divergence (MACD)**
- **Exponential Moving Average (EMA)**

====================================================================================
### File Dependencies:
- **CSV Data File:** `sl_training_Close.csv`
  - Format: `Date, Close Price`
====================================================================================
### Implementation Details:

- **Date Class:**
  - Parses and formats date strings from CSV files.

- **TimeSeries Struct:**
  - Stores vectors of dates and corresponding stock prices.

- **Core Functions:**
  - `readCSVFile()`: Reads stock data from CSV.
  - `findRSI()`: Computes the RSI for a given period.
  - `findEMA()`: Computes the Exponential Moving Average.
  - `findMACD()`: Calculates MACD from short/long EMAs.
  - `findSignalLine()`: Computes the MACD signal line.
  - `findMACDHistogram()`: Derives the MACD histogram.
  - `findSVI()`: (Placeholder) Computes Stochastic Volatility Index.
  - `getSignal()`: Determines buy/sell/hold recommendation.

- **Main Execution (`main()`):**
  - Defines periods for RSI, EMA, and MACD calculations.
  - Calls functions to compute indicators.
  - Prints results in a structured format.
*/
/*
Global Variables in `main()`
 - `const std::string dataDirectory = "sl_training_Close.csv";`
   - Stores the filename of the dataset to be processed.

 - `const int period = 5;`
   - The period length used for RSI and SVI calculations.

 - `const int shortTermPeriod = 12;`
   - The short EMA period for MACD calculation.

 - `const int longTermPeriod = 26;`
   - The long EMA period for MACD calculation.

 - `const int signalPeriod = 9;`
   - The period for the MACD signal line.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>
#include <sstream>
#include <memory>
#include <filesystem>
#include <algorithm>

namespace fs = std::filesystem;

class Date {
public:
    int year, month, day;

    Date(const std::string& dateStr) {
        std::stringstream ss(dateStr);
        char dash;
        ss >> year >> dash >> month >> dash >> day;
        if (ss.fail()) {
            throw std::invalid_argument("Invalid date format " + dateStr);
        }
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << std::setfill('0')
            << std::setw(2) << year << "-"
            << std::setw(2) << month << "-"
            << std::setw(2) << day;
        return oss.str();
    }
};

struct TimeSeries {
    std::vector<Date> dates;
    std::vector<double> values;
};

double stringToDouble(const std::string& str) {
    std::string cleanedStr = str;
    cleanedStr.erase(std::remove(cleanedStr.begin(), cleanedStr.end(), ','), cleanedStr.end());
    try {
        return std::stod(cleanedStr);
    } catch (const std::invalid_argument& e) {
        throw std::invalid_argument("Error converting to double " + str);
    }
}

TimeSeries readCSVFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Unable to open file " + filename);
    }

    TimeSeries ts;
    std::string line;
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string dateStr, valueStr;

        if (std::getline(ss, dateStr, ',') && std::getline(ss, valueStr)) {
            try {
                ts.dates.push_back(Date(dateStr));
                ts.values.push_back(stringToDouble(valueStr));
            } catch (const std::exception& e) {
                std::cerr << "Error for line " << line << "\n" << e.what() << std::endl;
                continue;
            }
        }
    }
    return ts;
}

std::vector<double> findRSI(const TimeSeries& data, int period) {
    std::vector<double> rsi(data.values.size(), std::numeric_limits<double>::quiet_NaN());

    double gain = 0.0;
    double loss = 0.0;

    for (int i = 1; i <= period; i++) {
        double change = data.values[i] - data.values[i - 1];
        if (change > 0) gain += change;
        else loss += -change;
    }

    double avgGain = gain / period;
    double avgLoss = loss / period;

    for (size_t i = period; i < data.values.size(); ++i) {
        double change = data.values[i] - data.values[i - 1];
        double currentGain = (change > 0) ? change : 0.0;
        double currentLoss = (change < 0) ? -change : 0.0;

        avgGain = ((avgGain * (period - 1)) + currentGain) / period;
        avgLoss = ((avgLoss * (period - 1)) + currentLoss) / period;

        double rs = (avgLoss == 0) ? 100 : avgGain / avgLoss;
        rsi[i] = 100.0 - (100.0 / (1.0 + rs));
    }

    return rsi;
}

std::vector<double> findEMA(const std::vector<double>& data, int period) {
    std::vector<double> ema(data.size(), std::numeric_limits<double>::quiet_NaN());

    double sum = 0.0;
    for (size_t i = 0; i < period; ++i) {
        sum += data[i];
    }

    double SMA = sum / period;

    for (size_t i = 0; i < period; ++i) {
        ema[i] = SMA;
    }

    double alpha = 2.0 / (1.0 + period);
    for (size_t i = period; i < data.size(); ++i) {
        ema[i] = data[i] * alpha + ema[i - 1] * (1 - alpha);
    }

    return ema;
}

std::vector<double> findMACD(const TimeSeries& data, int shortTermPeriod, int longTermPeriod) {
    std::vector<double> shortTermEMA = findEMA(data.values, shortTermPeriod);
    std::vector<double> longTermEMA = findEMA(data.values, longTermPeriod);

    std::vector<double> macd(data.values.size(), std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < data.values.size(); ++i) {
        if (std::isnan(shortTermEMA[i]) || std::isnan(longTermEMA[i])) {
            continue;
        }
        macd[i] = shortTermEMA[i] - longTermEMA[i];
    }

    return macd;
}

std::vector<double> findSignalLine(const std::vector<double>& macd, int signalPeriod) {
    return findEMA(macd, signalPeriod);
}

std::vector<double> findMACDHistogram(const std::vector<double>& macd, const std::vector<double>& signalLine) {
    std::vector<double> histogram(macd.size(), std::numeric_limits<double>::quiet_NaN());

    for (size_t i = 0; i < macd.size(); ++i) {
        if (std::isnan(macd[i]) || std::isnan(signalLine[i])) {
            continue;
        }
        histogram[i] = macd[i] - signalLine[i];
    }

    return histogram;
}

std::string getSignal(double rsi, double svi, double histogram) {
    if (std::isnan(rsi) || std::isnan(svi) || std::isnan(histogram)) {
        return "HOLD";
    }

    if (rsi < 30 && svi < 1.0 && histogram > 0) return "BUY";
    else if (rsi > 70 && svi > 1.0 && histogram <= 0) return "SELL";
    else return "HOLD";
}

std::vector<double> findSVI(const TimeSeries& data, int period) {
    std::vector<double> svi(data.values.size(), 1.0);
    return svi;
}

void printRSIandSVI(const TimeSeries& ts, const std::vector<double>& rsi, const std::vector<double>& svi, const std::vector<double>& macd, const std::vector<double>& signalLine, const std::vector<double>& histogram) {
    std::cout << std::setw(12) << "Date"
              << std::setw(12) << "RSI"
              << std::setw(12) << "SVI"
              << std::setw(12) << "MACD"
              << std::setw(12) << "Signal Line"
              << std::setw(12) << "Histogram"
              << std::setw(10) << "Signal" << std::endl;
    std::cout << std::string(70, '-') << std::endl;

    for (size_t i = 0; i < ts.dates.size(); ++i) {
        std::cout << std::setw(12) << ts.dates[i].toString();

        if (std::isnan(rsi[i])) std::cout << std::setw(12) << "NaN";
        else std::cout << std::setw(12) << std::fixed << std::setprecision(2) << rsi[i];

        if (std::isnan(svi[i])) std::cout << std::setw(12) << "NaN";
        else std::cout << std::setw(12) << std::fixed << std::setprecision(4) << svi[i];

        if (std::isnan(macd[i])) std::cout << std::setw(12) << "NaN";
        else std::cout << std::setw(12) << std::fixed << std::setprecision(4) << macd[i];

        if (std::isnan(signalLine[i])) std::cout << std::setw(12) << "NaN";
        else std::cout << std::setw(12) << std::fixed << std::setprecision(4) << signalLine[i];

        if (std::isnan(histogram[i])) std::cout << std::setw(12) << "NaN";
        else std::cout << std::setw(12) << std::fixed << std::setprecision(4) << histogram[i];

        std::cout << std::setw(10) << getSignal(rsi[i], svi[i], histogram[i]) << std::endl;
    }
}

int main() {
    const std::string dataDirectory = "sl_training_Close.csv";
    const int period = 5;
    const int shortTermPeriod = 12;
    const int longTermPeriod = 26;
    const int signalPeriod = 9;

    TimeSeries ts = readCSVFile(dataDirectory);
    std::vector<double> rsi = findRSI(ts, period);
    std::vector<double> svi = findSVI(ts, period);

    std::vector<double> macd = findMACD(ts, shortTermPeriod, longTermPeriod);
    std::vector<double> signalLine = findSignalLine(macd, signalPeriod);
    std::vector<double> histogram = findMACDHistogram(macd, signalLine);

    printRSIandSVI(ts, rsi, svi, macd, signalLine, histogram);

}
