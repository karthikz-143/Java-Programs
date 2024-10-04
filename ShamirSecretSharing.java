import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.JsonElement;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class ShamirSecretSharing {

    // Decoding y-value from different bases
    public static int decodeValue(String baseStr, String valueStr) {
        int base = Integer.parseInt(baseStr);
        return Integer.parseInt(valueStr, base);
    }

    // Extract (x, y) points from JSON
    public static void extractXYValues(JsonObject json, ArrayList<Double> xList, ArrayList<Double> yList) {
        for (Map.Entry<String, JsonElement> entry : json.entrySet()) {
            String key = entry.getKey();
            if (key.equals("keys")) {
                continue;
            }

            JsonObject point = entry.getValue().getAsJsonObject();
            int x = Integer.parseInt(key);
            int y = decodeValue(point.get("base").getAsString(), point.get("value").getAsString());

            xList.add((double) x);
            yList.add((double) y);
        }
    }

    // Generate combinations of indices to find subsets of k points
    public static List<int[]> combinations(int n, int k) {
        List<int[]> result = new ArrayList<>();
        combinationsHelper(result, new int[k], 0, n - 1, 0);
        return result;
    }

    // Helper method for generating combinations
    public static void combinationsHelper(List<int[]> result, int[] data, int start, int end, int index) {
        if (index == data.length) {
            result.add(data.clone());
            return;
        }
        for (int i = start; i <= end; i++) {
            data[index] = i;
            combinationsHelper(result, data, i + 1, end, index + 1);
        }
    }

    // Lagrange interpolation to find the polynomial coefficients
    public static double[] interpolatePolynomial(double[] x, double[] y) {
        int n = x.length;
        double[] coefficients = new double[n];
        for (int i = 0; i < n; i++) {
            double term = y[i];
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    term /= (x[i] - x[j]);
                }
            }
            for (int j = 0; j < n; j++) {
                double product = 1.0;
                for (int l = 0; l < n; l++) {
                    if (l != i && l != j) {
                        product *= (-x[l]);
                    }
                }
                coefficients[j] += term * product;
            }
        }
        return coefficients;
    }

    // Check if an index is part of the current subset
    public static boolean isInSubset(int index, int[] subset) {
        for (int i : subset) {
            if (i == index) {
                return true;
            }
        }
        return false;
    }

    // Evaluate the polynomial at a given x value
    public static double evaluatePolynomial(double[] polynomial, double x) {
        double result = 0;
        for (int i = polynomial.length - 1; i >= 0; i--) {
            result = result * x + polynomial[i];
        }
        return result;
    }

    // Find wrong points that don't lie on the curve
    public static List<int[]> findWrongPoints(double[] xArray, double[] yArray, int k) {
        List<int[]> wrongPoints = new ArrayList<>();
        for (int[] subset : combinations(xArray.length, k)) {
            double[] xSubset = new double[k];
            double[] ySubset = new double[k];
            for (int i = 0; i < k; i++) {
                xSubset[i] = xArray[subset[i]];
                ySubset[i] = yArray[subset[i]];
            }
            double[] polynomial = interpolatePolynomial(xSubset, ySubset);
            for (int i = 0; i < xArray.length; i++) {
                if (!isInSubset(i, subset)) {
                    double expectedY = evaluatePolynomial(polynomial, xArray[i]);
                    if (Math.abs(expectedY - yArray[i]) > 1e-6) {
                        wrongPoints.add(new int[]{(int) xArray[i], (int) yArray[i]});
                    }
                }
            }
            if (!wrongPoints.isEmpty()) {
                return wrongPoints;
            }
        }
        return wrongPoints;
    }

    // Compute the constant term 'c' using Lagrange interpolation
    public static double qlConstantTerm(double[] rs, double[] ry) {
        int n = rs.length;
        double x = 0.0;
        double[][] mlocx = new double[n][n];
        double[][] mloci = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                mlocx[i][j] = rs[i];
                mloci[i][j] = rs[i];
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    mlocx[i][j] = 1.0;
                } else {
                    mlocx[i][j] = -mlocx[i][j] + x;
                }
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    mloci[i][j] = 1.0;
                } else {
                    mloci[i][j] = -mloci[i][j] + rs[j];
                }
            }
        }
        double[] px = new double[n];
        double[] pi = new double[n];
        for (int i = 0; i < n; i++) {
            px[i] = 1.0;
            pi[i] = 1.0;
            for (int j = 0; j < n; j++) {
                px[i] *= mlocx[i][j];
                pi[i] *= mloci[i][j];
            }
        }
        double[] polyvect = new double[n];
        for (int i = 0; i < n; i++) {
            polyvect[i] = px[i] / pi[i];
        }
        double qlloc = 0.0;
        for (int i = 0; i < n; i++) {
            qlloc += ry[i] * polyvect[i];
        }
        return qlloc;
    }

    // Main function to parse the JSON, find wrong points, and calculate the constant term
    public static void main(String[] args) {
        String jsonFilePath = "input.json";
        try (FileReader reader = new FileReader(jsonFilePath)) {
            JsonObject json = JsonParser.parseReader(reader).getAsJsonObject();
            ArrayList<Double> xList = new ArrayList<>();
            ArrayList<Double> yList = new ArrayList<>();
            extractXYValues(json, xList, yList);
            double[] xArray = xList.stream().mapToDouble(Double::doubleValue).toArray();
            double[] yArray = yList.stream().mapToDouble(Double::doubleValue).toArray();
            int k = json.get("keys").getAsJsonObject().get("k").getAsInt();
            List<int[]> wrongPoints = findWrongPoints(xArray, yArray, k);
            if (wrongPoints.isEmpty()) {
                System.out.println("No wrong points found.");
            } else {
                for (int[] point : wrongPoints) {
                    System.out.println("Wrong point found: x = " + point[0] + ", y = " + point[1]);
                }
            }
            double constantTerm = qlConstantTerm(xArray, yArray);
            System.out.println("Constant Term: " + constantTerm);
        } catch (IOException e) {
            System.out.println("Error reading the JSON file: " + e.getMessage());
        }
    }
}
