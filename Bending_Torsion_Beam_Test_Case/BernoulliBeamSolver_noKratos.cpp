// Standalone 3D Euler-Bernoulli beam solver, no Kratos headers required.
//
// Example:
//   g++ -std=c++17 -O2 BernoulliBeamSolver_noKratos.cpp -o BernoulliBeamSolver_noKratos
//   ./BernoulliBeamSolver_noKratos --project Bending_Torsion_Beam_Test_Case/beam_geometry/ProjectParameters.json
//
// Useful overrides:
//   --mdpa path/to/model.mdpa
//   --materials path/to/StructuralMaterials.json
//   --output BernoulliBeamSolverOutput_vtk/BeamModelPart_0_0.vtk
//   --load node_id fx fy fz mx my mz
//   --fix-node node_id
//   --prescribe node_id dof value       dof = ux|uy|uz|rx|ry|rz

#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

namespace {

constexpr std::size_t DofsPerNode = 6;
constexpr double Pi = 3.141592653589793238462643383279502884;

struct Json {
    using Array = std::vector<Json>;
    using Object = std::map<std::string, Json>;
    std::variant<std::nullptr_t, bool, double, std::string, Array, Object> value = nullptr;

    bool IsNull() const { return std::holds_alternative<std::nullptr_t>(value); }
    bool IsBool() const { return std::holds_alternative<bool>(value); }
    bool IsNumber() const { return std::holds_alternative<double>(value); }
    bool IsString() const { return std::holds_alternative<std::string>(value); }
    bool IsArray() const { return std::holds_alternative<Array>(value); }
    bool IsObject() const { return std::holds_alternative<Object>(value); }

    bool Has(const std::string& key) const
    {
        return IsObject() && std::get<Object>(value).count(key) > 0;
    }

    const Json& operator[](const std::string& key) const
    {
        static const Json empty;
        if (!IsObject()) return empty;
        const auto& object = std::get<Object>(value);
        const auto it = object.find(key);
        return it == object.end() ? empty : it->second;
    }

    const Json& operator[](const std::size_t index) const
    {
        static const Json empty;
        if (!IsArray()) return empty;
        const auto& array = std::get<Array>(value);
        return index < array.size() ? array[index] : empty;
    }

    std::size_t Size() const
    {
        if (IsArray()) return std::get<Array>(value).size();
        if (IsObject()) return std::get<Object>(value).size();
        return 0;
    }

    bool AsBool(const bool fallback = false) const
    {
        return IsBool() ? std::get<bool>(value) : fallback;
    }

    double AsDouble(const double fallback = 0.0) const
    {
        return IsNumber() ? std::get<double>(value) : fallback;
    }

    std::string AsString(const std::string& fallback = "") const
    {
        return IsString() ? std::get<std::string>(value) : fallback;
    }

    const Array& AsArray() const
    {
        static const Array empty;
        return IsArray() ? std::get<Array>(value) : empty;
    }
};

std::string ReadFile(const std::filesystem::path& file_name)
{
    std::ifstream file(file_name);
    if (!file) {
        throw std::runtime_error("Could not open file: " + file_name.string());
    }
    return std::string((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
}

std::string StripJsonComments(const std::string& text)
{
    std::string out;
    bool in_string = false;
    bool escaped = false;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char c = text[i];
        if (in_string) {
            out.push_back(c);
            escaped = (!escaped && c == '\\');
            if (!escaped && c == '"') in_string = false;
            if (c != '\\') escaped = false;
        } else if (c == '"') {
            in_string = true;
            out.push_back(c);
        } else if (c == '/' && i + 1 < text.size() && text[i + 1] == '/') {
            while (i < text.size() && text[i] != '\n') ++i;
            out.push_back('\n');
        } else {
            out.push_back(c);
        }
    }
    return out;
}

class JsonParser {
public:
    explicit JsonParser(std::string text) : mText(std::move(text)) {}

    Json Parse()
    {
        SkipWhitespace();
        Json result = ParseValue();
        SkipWhitespace();
        return result;
    }

private:
    Json ParseValue()
    {
        SkipWhitespace();
        if (Match('{')) return ParseObject();
        if (Match('[')) return ParseArray();
        if (Peek() == '"') return Json{ParseString()};
        if (StartsWith("true")) {
            mPosition += 4;
            return Json{true};
        }
        if (StartsWith("false")) {
            mPosition += 5;
            return Json{false};
        }
        if (StartsWith("null")) {
            mPosition += 4;
            return Json{nullptr};
        }
        return Json{ParseNumber()};
    }

    Json ParseObject()
    {
        Json::Object object;
        SkipWhitespace();
        if (Match('}')) return Json{object};
        while (true) {
            SkipWhitespace();
            const std::string key = ParseString();
            Require(':');
            object[key] = ParseValue();
            SkipWhitespace();
            if (Match('}')) break;
            Require(',');
        }
        return Json{object};
    }

    Json ParseArray()
    {
        Json::Array array;
        SkipWhitespace();
        if (Match(']')) return Json{array};
        while (true) {
            array.push_back(ParseValue());
            SkipWhitespace();
            if (Match(']')) break;
            Require(',');
        }
        return Json{array};
    }

    std::string ParseString()
    {
        Require('"');
        std::string result;
        while (mPosition < mText.size()) {
            const char c = mText[mPosition++];
            if (c == '"') break;
            if (c == '\\' && mPosition < mText.size()) {
                const char e = mText[mPosition++];
                if (e == 'n') result.push_back('\n');
                else if (e == 't') result.push_back('\t');
                else result.push_back(e);
            } else {
                result.push_back(c);
            }
        }
        return result;
    }

    double ParseNumber()
    {
        const std::size_t start = mPosition;
        while (mPosition < mText.size()) {
            const char c = mText[mPosition];
            if (!std::isdigit(c) && c != '-' && c != '+' && c != '.' && c != 'e' && c != 'E') break;
            ++mPosition;
        }
        if (start == mPosition) {
            throw std::runtime_error("Invalid JSON near: " + mText.substr(mPosition, 40));
        }
        return std::stod(mText.substr(start, mPosition - start));
    }

    void SkipWhitespace()
    {
        while (mPosition < mText.size() && std::isspace(static_cast<unsigned char>(mText[mPosition]))) ++mPosition;
    }

    bool StartsWith(const std::string& word) const
    {
        return mText.compare(mPosition, word.size(), word) == 0;
    }

    char Peek() const
    {
        return mPosition < mText.size() ? mText[mPosition] : '\0';
    }

    bool Match(const char c)
    {
        SkipWhitespace();
        if (Peek() == c) {
            ++mPosition;
            return true;
        }
        return false;
    }

    void Require(const char c)
    {
        if (!Match(c)) {
            throw std::runtime_error(std::string("Expected '") + c + "' in JSON.");
        }
    }

    std::string mText;
    std::size_t mPosition = 0;
};

Json ReadJson(const std::filesystem::path& file_name)
{
    return JsonParser(StripJsonComments(ReadFile(file_name))).Parse();
}

class Expression {
public:
    Expression(std::string expression, const double x, const double y, const double z, const double t)
        : mExpression(std::move(expression)), mX(x), mY(y), mZ(z), mT(t)
    {
    }

    double Evaluate()
    {
        const double value = ParseExpression();
        SkipWhitespace();
        return value;
    }

private:
    double ParseExpression()
    {
        double value = ParseTerm();
        while (true) {
            SkipWhitespace();
            if (Match('+')) value += ParseTerm();
            else if (Match('-')) value -= ParseTerm();
            else return value;
        }
    }

    double ParseTerm()
    {
        double value = ParsePower();
        while (true) {
            SkipWhitespace();
            if (Match('*')) value *= ParsePower();
            else if (Match('/')) value /= ParsePower();
            else return value;
        }
    }

    double ParsePower()
    {
        double value = ParseUnary();
        SkipWhitespace();
        if (Match('^')) value = std::pow(value, ParsePower());
        return value;
    }

    double ParseUnary()
    {
        SkipWhitespace();
        if (Match('+')) return ParseUnary();
        if (Match('-')) return -ParseUnary();
        return ParsePrimary();
    }

    double ParsePrimary()
    {
        SkipWhitespace();
        if (Match('(')) {
            const double value = ParseExpression();
            Require(')');
            return value;
        }
        if (std::isdigit(Peek()) || Peek() == '.') return ParseNumber();
        if (std::isalpha(static_cast<unsigned char>(Peek())) || Peek() == '_') {
            const std::string name = ParseIdentifier();
            SkipWhitespace();
            if (Match('(')) {
                const double argument = ParseExpression();
                Require(')');
                if (name == "sin") return std::sin(argument);
                if (name == "cos") return std::cos(argument);
                if (name == "tan") return std::tan(argument);
                if (name == "sqrt") return std::sqrt(argument);
                if (name == "exp") return std::exp(argument);
                if (name == "abs") return std::abs(argument);
                throw std::runtime_error("Unsupported function in expression: " + name);
            }
            if (name == "x" || name == "X") return mX;
            if (name == "y" || name == "Y") return mY;
            if (name == "z" || name == "Z") return mZ;
            if (name == "t") return mT;
            if (name == "pi") return Pi;
            throw std::runtime_error("Unsupported symbol in expression: " + name);
        }
        throw std::runtime_error("Invalid expression near: " + mExpression.substr(mPosition, 30));
    }

    double ParseNumber()
    {
        const std::size_t start = mPosition;
        while (mPosition < mExpression.size()) {
            const char c = mExpression[mPosition];
            if (!std::isdigit(c) && c != '.' && c != 'e' && c != 'E' && c != '+' && c != '-') break;
            if ((c == '+' || c == '-') && mPosition != start) {
                const char previous = mExpression[mPosition - 1];
                if (previous != 'e' && previous != 'E') break;
            }
            ++mPosition;
        }
        return std::stod(mExpression.substr(start, mPosition - start));
    }

    std::string ParseIdentifier()
    {
        const std::size_t start = mPosition;
        while (mPosition < mExpression.size()) {
            const char c = mExpression[mPosition];
            if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_') break;
            ++mPosition;
        }
        return mExpression.substr(start, mPosition - start);
    }

    void SkipWhitespace()
    {
        while (mPosition < mExpression.size() && std::isspace(static_cast<unsigned char>(mExpression[mPosition]))) ++mPosition;
    }

    char Peek() const
    {
        return mPosition < mExpression.size() ? mExpression[mPosition] : '\0';
    }

    bool Match(const char c)
    {
        SkipWhitespace();
        if (Peek() == c) {
            ++mPosition;
            return true;
        }
        return false;
    }

    void Require(const char c)
    {
        if (!Match(c)) throw std::runtime_error(std::string("Expected '") + c + "' in expression.");
    }

    std::string mExpression;
    std::size_t mPosition = 0;
    double mX = 0.0, mY = 0.0, mZ = 0.0, mT = 0.0;
};

struct Node {
    int id = 0;
    std::array<double, 3> x{};
    std::array<double, 6> values{};
    std::array<double, 6> loads{};
    std::array<bool, 6> fixed{};
};

struct Element {
    int id = 0;
    int node1 = 0;
    int node2 = 0;
};

struct Condition {
    int id = 0;
    std::vector<int> node_ids;
};

struct Model {
    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::vector<Condition> conditions;
    std::unordered_map<int, std::size_t> node_index_by_id;
    std::unordered_map<int, std::size_t> condition_index_by_id;
    std::unordered_map<std::string, std::vector<int>> sub_model_part_nodes;
    std::unordered_map<std::string, std::vector<int>> sub_model_part_conditions;
    std::unordered_map<std::string, double> material_values;
};

struct Material {
    double E = 206.9e9;
    double nu = 0.29;
    double A = 1.0;
    double Iy = 1.0;
    double Iz = 1.0;
    double J = 1.0;
    double G() const { return E / (2.0 * (1.0 + nu)); }
};

struct Settings {
    std::filesystem::path project_file;
    std::filesystem::path mdpa_file;
    std::filesystem::path materials_file;
    std::filesystem::path output_file = "BernoulliBeamSolverOutput_vtk/noKratos_Beam.vtk";
    Json project;
};

using Matrix12 = std::array<std::array<double, 12>, 12>;
using SparseRows = std::vector<std::unordered_map<int, double>>;

std::string StripMdpaComment(const std::string& line)
{
    const auto pos = line.find("//");
    return pos == std::string::npos ? line : line.substr(0, pos);
}

std::string Trim(const std::string& text)
{
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

std::vector<std::string> SplitWords(const std::string& line)
{
    std::istringstream stream(line);
    std::vector<std::string> words;
    std::string word;
    while (stream >> word) words.push_back(word);
    return words;
}

bool ExistsWithOptionalMdpa(const std::filesystem::path& path)
{
    std::filesystem::path mdpa_path = path;
    mdpa_path.replace_extension(".mdpa");
    return std::filesystem::exists(path) || std::filesystem::exists(mdpa_path);
}

std::filesystem::path ResolveInputPath(
    const std::filesystem::path& project_directory,
    const std::filesystem::path& relative_path)
{
    if (relative_path.is_absolute()) return relative_path;
    const auto from_project_dir = project_directory / relative_path;
    if (ExistsWithOptionalMdpa(from_project_dir)) return from_project_dir;
    const auto from_case_dir = project_directory.parent_path() / relative_path;
    if (ExistsWithOptionalMdpa(from_case_dir)) return from_case_dir;
    return from_project_dir;
}

std::filesystem::path WithoutMdpaExtension(std::filesystem::path path)
{
    if (path.extension() == ".mdpa") path.replace_extension("");
    return path;
}

Settings ReadSettings(const std::filesystem::path& project_file)
{
    Settings settings;
    settings.project_file = project_file;
    settings.project = ReadJson(project_file);

    const Json& solver = settings.project["solver_settings"];
    if (!solver.Has("model_import_settings") && solver.Has("solvers") && solver["solvers"].Has("beam_structure")) {
        std::filesystem::path nested = solver["solvers"]["beam_structure"]["solver_wrapper_settings"]["input_file"].AsString();
        if (nested.extension().empty()) nested.replace_extension(".json");
        return ReadSettings(ResolveInputPath(project_file.parent_path(), nested));
    }

    settings.mdpa_file = ResolveInputPath(
        project_file.parent_path(),
        solver["model_import_settings"]["input_filename"].AsString());
    settings.materials_file = ResolveInputPath(
        project_file.parent_path(),
        solver["material_import_settings"]["materials_filename"].AsString());

    const Json& vtk_outputs = settings.project["output_processes"]["vtk_output"];
    if (vtk_outputs.IsArray() && vtk_outputs.Size() > 0 && vtk_outputs[0]["Parameters"].Has("output_path")) {
        std::filesystem::path output_path = vtk_outputs[0]["Parameters"]["output_path"].AsString();
        if (output_path.extension() != ".vtk") output_path /= "BeamModelPart.vtk";
        settings.output_file = ResolveInputPath(project_file.parent_path().parent_path(), output_path);
    }
    return settings;
}

Model ReadMdpa(const std::filesystem::path& mdpa_file)
{
    std::ifstream file(WithoutMdpaExtension(mdpa_file).replace_extension(".mdpa"));
    if (!file) {
        file.open(mdpa_file);
    }
    if (!file) throw std::runtime_error("Could not open mdpa file: " + mdpa_file.string());

    Model model;
    std::string line;
    std::string block;
    std::string current_sub_model_part;
    bool in_sub_nodes = false;
    bool in_sub_conditions = false;

    while (std::getline(file, line)) {
        const auto words = SplitWords(Trim(StripMdpaComment(line)));
        if (words.empty()) continue;

        if (words[0] == "Begin") {
            if (words.size() >= 2 && words[1] == "Nodes") block = "Nodes";
            else if (words.size() >= 2 && words[1] == "Elements") block = "Elements";
            else if (words.size() >= 2 && words[1] == "Conditions") block = "Conditions";
            else if (words.size() >= 2 && words[1] == "Properties") block = "Properties";
            else if (words.size() >= 3 && words[1] == "SubModelPart") {
                current_sub_model_part = words[2];
                model.sub_model_part_nodes[current_sub_model_part];
                model.sub_model_part_conditions[current_sub_model_part];
                block = "SubModelPart";
            } else if (words.size() >= 2 && words[1] == "SubModelPartNodes") in_sub_nodes = true;
            else if (words.size() >= 2 && words[1] == "SubModelPartConditions") in_sub_conditions = true;
            continue;
        }

        if (words[0] == "End") {
            if (words.size() >= 2 && words[1] == "Nodes") block.clear();
            else if (words.size() >= 2 && words[1] == "Elements") block.clear();
            else if (words.size() >= 2 && words[1] == "Conditions") block.clear();
            else if (words.size() >= 2 && words[1] == "Properties") block.clear();
            else if (words.size() >= 2 && words[1] == "SubModelPartNodes") in_sub_nodes = false;
            else if (words.size() >= 2 && words[1] == "SubModelPartConditions") in_sub_conditions = false;
            else if (words.size() >= 2 && words[1] == "SubModelPart") {
                current_sub_model_part.clear();
                block.clear();
            }
            continue;
        }

        if (block == "Nodes" && words.size() >= 4) {
            Node node;
            node.id = std::stoi(words[0]);
            node.x = {std::stod(words[1]), std::stod(words[2]), std::stod(words[3])};
            model.node_index_by_id[node.id] = model.nodes.size();
            model.nodes.push_back(node);
        } else if (block == "Elements" && words.size() >= 4) {
            model.elements.push_back({std::stoi(words[0]), std::stoi(words[2]), std::stoi(words[3])});
        } else if (block == "Conditions" && words.size() >= 3) {
            Condition condition;
            condition.id = std::stoi(words[0]);
            for (std::size_t i = 2; i < words.size(); ++i) condition.node_ids.push_back(std::stoi(words[i]));
            model.condition_index_by_id[condition.id] = model.conditions.size();
            model.conditions.push_back(condition);
        } else if (block == "Properties" && words.size() >= 2) {
            try {
                model.material_values[words[0]] = std::stod(words[1]);
            } catch (const std::exception&) {
                // Non-scalar mdpa property entries are intentionally ignored by this minimal reader.
            }
        } else if (!current_sub_model_part.empty() && in_sub_nodes) {
            model.sub_model_part_nodes[current_sub_model_part].push_back(std::stoi(words[0]));
        } else if (!current_sub_model_part.empty() && in_sub_conditions) {
            model.sub_model_part_conditions[current_sub_model_part].push_back(std::stoi(words[0]));
        }
    }
    if (model.nodes.empty() || model.elements.empty()) {
        throw std::runtime_error("The mdpa must contain nodes and 2-node beam elements.");
    }
    return model;
}

void ApplyMdpaMaterialOverrides(const Model& model, Material& material)
{
    auto assign_if_found = [&model](const std::string& key, double& rValue) {
        const auto it = model.material_values.find(key);
        if (it != model.material_values.end()) rValue = it->second;
    };
    assign_if_found("YOUNG_MODULUS", material.E);
    assign_if_found("POISSON_RATIO", material.nu);
    assign_if_found("CROSS_AREA", material.A);
    assign_if_found("I22", material.Iy);
    assign_if_found("I33", material.Iz);
    assign_if_found("TORSIONAL_INERTIA", material.J);
}

Material ReadMaterial(const std::filesystem::path& material_file, const Model& model)
{
    Material material;
    if (!std::filesystem::exists(material_file)) {
        std::cerr << "Warning: material file not found; using default material.\n";
        ApplyMdpaMaterialOverrides(model, material);
        return material;
    }
    const Json json = ReadJson(material_file);
    const Json& variables = json["properties"][0]["Material"]["Variables"];
    material.E = variables["YOUNG_MODULUS"].AsDouble(material.E);
    material.nu = variables["POISSON_RATIO"].AsDouble(material.nu);
    material.A = variables["CROSS_AREA"].AsDouble(material.A);
    material.Iy = variables["I22"].AsDouble(material.Iy);
    material.Iz = variables["I33"].AsDouble(material.Iz);
    material.J = variables["TORSIONAL_INERTIA"].AsDouble(material.J);
    ApplyMdpaMaterialOverrides(model, material);
    return material;
}

std::vector<int> NodeIdsForModelPart(const Model& model, const std::string& model_part_name)
{
    const auto dot = model_part_name.find('.');
    const std::string short_name = dot == std::string::npos ? model_part_name : model_part_name.substr(model_part_name.rfind('.') + 1);
    const auto it = model.sub_model_part_nodes.find(short_name);
    if (model_part_name == "Structure" || model_part_name == "BeamModelPart" || it == model.sub_model_part_nodes.end()) {
        std::vector<int> ids;
        ids.reserve(model.nodes.size());
        for (const auto& node : model.nodes) ids.push_back(node.id);
        return ids;
    }
    return it->second;
}

std::vector<int> ConditionIdsForModelPart(const Model& model, const std::string& model_part_name)
{
    const std::string short_name = model_part_name.substr(model_part_name.rfind('.') + 1);
    const auto it = model.sub_model_part_conditions.find(short_name);
    if (it == model.sub_model_part_conditions.end()) {
        std::vector<int> ids;
        ids.reserve(model.conditions.size());
        for (const auto& condition : model.conditions) ids.push_back(condition.id);
        return ids;
    }
    return it->second;
}

int DofIndex(const std::string& dof)
{
    static const std::unordered_map<std::string, int> map = {
        {"ux", 0}, {"uy", 1}, {"uz", 2}, {"rx", 3}, {"ry", 4}, {"rz", 5},
        {"DISPLACEMENT_X", 0}, {"DISPLACEMENT_Y", 1}, {"DISPLACEMENT_Z", 2},
        {"ROTATION_X", 3}, {"ROTATION_Y", 4}, {"ROTATION_Z", 5}};
    const auto it = map.find(dof);
    if (it == map.end()) throw std::runtime_error("Unknown dof: " + dof);
    return it->second;
}

double JsonValueAtNode(const Json& value, const Node& node, const double time)
{
    if (value.IsNumber()) return value.AsDouble();
    if (value.IsString()) return Expression(value.AsString(), node.x[0], node.x[1], node.x[2], time).Evaluate();
    throw std::runtime_error("JSON process value must be a number or expression string.");
}

void SetVectorValue(Node& node, const std::string& variable, const std::size_t component, const double value)
{
    if (variable == "DISPLACEMENT") node.values[component] = value;
    else if (variable == "ROTATION") node.values[component + 3] = value;
    else if (variable == "POINT_LOAD") node.loads[component] += value;
    else if (variable == "POINT_MOMENT") node.loads[component + 3] += value;
    else throw std::runtime_error("Unsupported variable in process: " + variable);
}

void ApplyAssignVectorVariableProcess(Model& model, const Json& parameters, const double time, const bool fixes_components)
{
    const std::string model_part_name = parameters["model_part_name"].AsString("Structure");
    const std::string variable = parameters["variable_name"].AsString();
    const Json& values = parameters["value"];
    const Json& constrained = parameters["constrained"];

    for (const int node_id : NodeIdsForModelPart(model, model_part_name)) {
        Node& node = model.nodes.at(model.node_index_by_id.at(node_id));
        for (std::size_t component = 0; component < 3; ++component) {
            SetVectorValue(node, variable, component, JsonValueAtNode(values[component], node, time));
            if (fixes_components && constrained.IsArray() && constrained[component].AsBool()) {
                if (variable == "DISPLACEMENT") node.fixed[component] = true;
                else if (variable == "ROTATION") node.fixed[component + 3] = true;
            }
        }
    }
}

void ApplyVectorByDirectionToConditionProcess(Model& model, const Json& parameters)
{
    const std::string model_part_name = parameters["model_part_name"].AsString();
    const std::string variable = parameters["variable_name"].AsString("POINT_LOAD");
    const double modulus = parameters["modulus"].AsDouble();
    std::array<double, 3> value{};
    for (std::size_t i = 0; i < 3; ++i) value[i] = modulus * parameters["direction"][i].AsDouble();

    for (const int condition_id : ConditionIdsForModelPart(model, model_part_name)) {
        const Condition& condition = model.conditions.at(model.condition_index_by_id.at(condition_id));
        const double share = 1.0 / static_cast<double>(condition.node_ids.size());
        for (const int node_id : condition.node_ids) {
            Node& node = model.nodes.at(model.node_index_by_id.at(node_id));
            for (std::size_t component = 0; component < 3; ++component) {
                SetVectorValue(node, variable, component, share * value[component]);
            }
        }
    }
}

void ApplyJsonProcesses(Model& model, const Json& project)
{
    if (!project.Has("processes")) return;

    if (project["solver_settings"]["domain_size"].AsDouble(3.0) == 2.0) {
        for (auto& node : model.nodes) {
            node.fixed[2] = true;
            node.fixed[4] = true;
        }
    }

    const double process_time = project["problem_data"]["end_time"].AsDouble(0.0);
    for (const auto& process : project["processes"]["constraints_process_list"].AsArray()) {
        const std::string process_name = process["process_name"].AsString();
        if (process_name == "AssignVectorVariableProcess") {
            ApplyAssignVectorVariableProcess(model, process["Parameters"], process_time, true);
        }
    }
    for (const auto& process : project["processes"]["loads_process_list"].AsArray()) {
        const std::string process_name = process["process_name"].AsString();
        if (process_name == "AssignVectorByDirectionToConditionProcess") {
            ApplyVectorByDirectionToConditionProcess(model, process["Parameters"]);
        } else if (process_name == "AssignVectorVariableProcess") {
            ApplyAssignVectorVariableProcess(model, process["Parameters"], process_time, false);
        }
    }
}

double Dot(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

std::array<double, 3> Cross(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

double Norm(const std::array<double, 3>& v)
{
    return std::sqrt(Dot(v, v));
}

std::array<double, 3> Normalize(const std::array<double, 3>& v)
{
    const double norm = Norm(v);
    if (norm <= std::numeric_limits<double>::epsilon()) throw std::runtime_error("Zero vector normalization.");
    return {v[0] / norm, v[1] / norm, v[2] / norm};
}

Matrix12 LocalStiffness(const Material& material, const double L)
{
    Matrix12 k{};
    const double EA_L = material.E * material.A / L;
    k[0][0] = EA_L; k[0][6] = -EA_L; k[6][0] = -EA_L; k[6][6] = EA_L;
    const double GJ_L = material.G() * material.J / L;
    k[3][3] = GJ_L; k[3][9] = -GJ_L; k[9][3] = -GJ_L; k[9][9] = GJ_L;

    auto add = [&k](const std::array<int, 4>& dofs, const std::array<std::array<double, 4>, 4>& kb) {
        for (std::size_t i = 0; i < 4; ++i) for (std::size_t j = 0; j < 4; ++j) k[dofs[i]][dofs[j]] += kb[i][j];
    };

    double c1 = 12.0 * material.E * material.Iz / std::pow(L, 3);
    double c2 = 6.0 * material.E * material.Iz / std::pow(L, 2);
    double c3 = 4.0 * material.E * material.Iz / L;
    double c4 = 2.0 * material.E * material.Iz / L;
    add({1, 5, 7, 11}, {{{c1, c2, -c1, c2}, {c2, c3, -c2, c4}, {-c1, -c2, c1, -c2}, {c2, c4, -c2, c3}}});

    c1 = 12.0 * material.E * material.Iy / std::pow(L, 3);
    c2 = 6.0 * material.E * material.Iy / std::pow(L, 2);
    c3 = 4.0 * material.E * material.Iy / L;
    c4 = 2.0 * material.E * material.Iy / L;
    add({2, 4, 8, 10}, {{{c1, -c2, -c1, -c2}, {-c2, c3, c2, c4}, {-c1, c2, c1, c2}, {-c2, c4, c2, c3}}});
    return k;
}

Matrix12 TransformToGlobal(const Matrix12& local_k, const Node& node1, const Node& node2)
{
    const auto ex = Normalize({node2.x[0] - node1.x[0], node2.x[1] - node1.x[1], node2.x[2] - node1.x[2]});
    std::array<double, 3> reference = {0.0, 0.0, 1.0};
    if (std::abs(Dot(ex, reference)) > 0.95) reference = {0.0, 1.0, 0.0};
    const auto ey = Normalize(Cross(reference, ex));
    const auto ez = Cross(ex, ey);
    const std::array<std::array<double, 3>, 3> R = {ex, ey, ez};

    Matrix12 global_k{};
    for (std::size_t a = 0; a < 12; ++a) {
        const std::size_t a_node = a / DofsPerNode;
        const std::size_t a_component = a % 3;
        const bool a_rotation = (a % DofsPerNode) >= 3;
        for (std::size_t b = 0; b < 12; ++b) {
            const std::size_t b_node = b / DofsPerNode;
            const std::size_t b_component = b % 3;
            const bool b_rotation = (b % DofsPerNode) >= 3;
            double value = 0.0;
            for (std::size_t i = 0; i < 3; ++i) {
                const std::size_t local_a = a_node * DofsPerNode + (a_rotation ? 3 + i : i);
                for (std::size_t j = 0; j < 3; ++j) {
                    const std::size_t local_b = b_node * DofsPerNode + (b_rotation ? 3 + j : j);
                    value += R[i][a_component] * local_k[local_a][local_b] * R[j][b_component];
                }
            }
            global_k[a][b] = value;
        }
    }
    return global_k;
}

void Assemble(const Model& model, const Material& material, SparseRows& K)
{
    K.assign(model.nodes.size() * DofsPerNode, {});
    for (const auto& element : model.elements) {
        const auto i1 = model.node_index_by_id.at(element.node1);
        const auto i2 = model.node_index_by_id.at(element.node2);
        const Node& node1 = model.nodes[i1];
        const Node& node2 = model.nodes[i2];
        const double L = Norm({node2.x[0] - node1.x[0], node2.x[1] - node1.x[1], node2.x[2] - node1.x[2]});
        const Matrix12 kg = TransformToGlobal(LocalStiffness(material, L), node1, node2);
        std::array<int, 12> dofs{};
        for (std::size_t d = 0; d < DofsPerNode; ++d) {
            dofs[d] = static_cast<int>(i1 * DofsPerNode + d);
            dofs[DofsPerNode + d] = static_cast<int>(i2 * DofsPerNode + d);
        }
        for (std::size_t i = 0; i < 12; ++i) {
            for (std::size_t j = 0; j < 12; ++j) K[dofs[i]][dofs[j]] += kg[i][j];
        }
    }
}

std::vector<double> SolveReducedBandedCholesky(const SparseRows& K, const Model& model)
{
    const std::size_t full_size = K.size();
    std::vector<int> full_to_free(full_size, -1);
    std::vector<int> free_to_full;
    std::vector<double> fixed_values(full_size, 0.0);
    std::vector<double> full_loads(full_size, 0.0);

    for (std::size_t n = 0; n < model.nodes.size(); ++n) {
        for (std::size_t d = 0; d < DofsPerNode; ++d) {
            const std::size_t id = n * DofsPerNode + d;
            fixed_values[id] = model.nodes[n].values[d];
            full_loads[id] = model.nodes[n].loads[d];
            if (!model.nodes[n].fixed[d]) {
                full_to_free[id] = static_cast<int>(free_to_full.size());
                free_to_full.push_back(static_cast<int>(id));
            }
        }
    }

    if (free_to_full.empty()) return fixed_values;

    std::size_t bandwidth = 0;
    for (std::size_t i_full = 0; i_full < full_size; ++i_full) {
        const int i_free = full_to_free[i_full];
        if (i_free < 0) continue;
        for (const auto& [j_full, value] : K[i_full]) {
            const int j_free = full_to_free[j_full];
            if (j_free >= 0 && value != 0.0) bandwidth = std::max(bandwidth, static_cast<std::size_t>(std::abs(i_free - j_free)));
        }
    }

    const std::size_t n = free_to_full.size();
    std::vector<double> A(n * (bandwidth + 1), 0.0);
    auto lower = [bandwidth, &A](std::size_t i, std::size_t j) -> double& { return A[i * (bandwidth + 1) + (i - j)]; };
    auto lower_value = [bandwidth, &A](std::size_t i, std::size_t j) -> double {
        return (i >= j && i - j <= bandwidth) ? A[i * (bandwidth + 1) + (i - j)] : 0.0;
    };

    std::vector<double> b(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const int full_i = free_to_full[i];
        b[i] = full_loads[full_i];
        for (const auto& [full_j, value] : K[full_i]) {
            const int j_free = full_to_free[full_j];
            if (j_free >= 0) lower(std::max<std::size_t>(i, j_free), std::min<std::size_t>(i, j_free)) = value;
            else b[i] -= value * fixed_values[full_j];
        }
    }

    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t j_begin = i > bandwidth ? i - bandwidth : 0;
        for (std::size_t j = j_begin; j <= i; ++j) {
            double sum = lower(i, j);
            const std::size_t k_begin = std::max(i > bandwidth ? i - bandwidth : 0, j > bandwidth ? j - bandwidth : 0);
            for (std::size_t k = k_begin; k < j; ++k) sum -= lower_value(i, k) * lower_value(j, k);
            if (i == j) {
                if (sum <= 0.0) throw std::runtime_error("Reduced stiffness is not positive definite. Check BCs.");
                lower(i, j) = std::sqrt(sum);
            } else {
                lower(i, j) = sum / lower(j, j);
            }
        }
    }

    std::vector<double> y(n, 0.0), x(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        double sum = b[i];
        const std::size_t j_begin = i > bandwidth ? i - bandwidth : 0;
        for (std::size_t j = j_begin; j < i; ++j) sum -= lower_value(i, j) * y[j];
        y[i] = sum / lower(i, i);
    }
    for (std::size_t reverse_i = 0; reverse_i < n; ++reverse_i) {
        const std::size_t i = n - 1 - reverse_i;
        double sum = y[i];
        const std::size_t j_end = std::min(n - 1, i + bandwidth);
        for (std::size_t j = i + 1; j <= j_end; ++j) sum -= lower_value(j, i) * x[j];
        x[i] = sum / lower(i, i);
    }

    std::vector<double> U = fixed_values;
    for (std::size_t i = 0; i < n; ++i) U[free_to_full[i]] = x[i];
    std::cout << "Solved reduced system with " << n << " free DOFs and half-bandwidth " << bandwidth << ".\n";
    return U;
}

void WriteVtk(const std::filesystem::path& output_file, const Model& model, const std::vector<double>& U)
{
    std::filesystem::create_directories(output_file.parent_path().empty() ? "." : output_file.parent_path());
    std::ofstream out(output_file);
    if (!out) throw std::runtime_error("Could not write VTK file: " + output_file.string());

    out << "# vtk DataFile Version 3.0\nBernoulli beam solution\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << model.nodes.size() << " double\n" << std::setprecision(12);
    for (const auto& node : model.nodes) {
        const std::size_t base = model.node_index_by_id.at(node.id) * DofsPerNode;
        out << node.x[0] + U[base] << ' ' << node.x[1] + U[base + 1] << ' ' << node.x[2] + U[base + 2] << '\n';
    }
    out << "CELLS " << model.elements.size() << ' ' << model.elements.size() * 3 << '\n';
    for (const auto& element : model.elements) {
        out << "2 " << model.node_index_by_id.at(element.node1) << ' ' << model.node_index_by_id.at(element.node2) << '\n';
    }
    out << "CELL_TYPES " << model.elements.size() << '\n';
    for (std::size_t i = 0; i < model.elements.size(); ++i) out << "3\n";
    out << "POINT_DATA " << model.nodes.size() << "\nVECTORS DISPLACEMENT double\n";
    for (std::size_t i = 0; i < model.nodes.size(); ++i) {
        const std::size_t base = i * DofsPerNode;
        out << U[base] << ' ' << U[base + 1] << ' ' << U[base + 2] << '\n';
    }
    out << "VECTORS ROTATION double\n";
    for (std::size_t i = 0; i < model.nodes.size(); ++i) {
        const std::size_t base = i * DofsPerNode;
        out << U[base + 3] << ' ' << U[base + 4] << ' ' << U[base + 5] << '\n';
    }
}

void PrintUsage(const char* executable)
{
    std::cout << "Usage: " << executable << " [ProjectParameters.json|mesh.mdpa] [options]\n"
              << "  --project file.json\n"
              << "  --mdpa file.mdpa\n"
              << "  --materials StructuralMaterials.json\n"
              << "  --output result.vtk\n"
              << "  --load node fx fy fz mx my mz\n"
              << "  --fix-node node\n"
              << "  --prescribe node ux|uy|uz|rx|ry|rz value\n";
}

} // namespace

int main(int argc, char** argv)
{
    try {
        Settings settings;
        std::filesystem::path project_file =
            "Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/ProjectParameters.json";
        std::filesystem::path mdpa_override;
        std::filesystem::path materials_override;
        std::filesystem::path output_override;
        std::vector<std::pair<int, std::array<double, 6>>> load_overrides;
        std::vector<int> fixed_node_overrides;
        std::vector<std::tuple<int, int, double>> prescribed_overrides;

        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--help" || arg == "-h") {
                PrintUsage(argv[0]);
                return 0;
            } else if (arg == "--project" && i + 1 < argc) {
                project_file = argv[++i];
            } else if (arg == "--mdpa" && i + 1 < argc) {
                mdpa_override = argv[++i];
            } else if (arg == "--materials" && i + 1 < argc) {
                materials_override = argv[++i];
            } else if (arg == "--output" && i + 1 < argc) {
                output_override = argv[++i];
            } else if (arg == "--load" && i + 7 < argc) {
                const int node_id = std::stoi(argv[++i]);
                std::array<double, 6> load{};
                for (double& v : load) v = std::stod(argv[++i]);
                load_overrides.emplace_back(node_id, load);
            } else if (arg == "--fix-node" && i + 1 < argc) {
                fixed_node_overrides.push_back(std::stoi(argv[++i]));
            } else if (arg == "--prescribe" && i + 3 < argc) {
                const int node_id = std::stoi(argv[++i]);
                const int dof = DofIndex(argv[++i]);
                const double value = std::stod(argv[++i]);
                prescribed_overrides.emplace_back(node_id, dof, value);
            } else if (std::filesystem::path(arg).extension() == ".json") {
                project_file = arg;
            } else {
                mdpa_override = arg;
            }
        }

        settings = ReadSettings(project_file);
        if (!mdpa_override.empty()) settings.mdpa_file = mdpa_override;
        if (!materials_override.empty()) settings.materials_file = materials_override;
        if (!output_override.empty()) settings.output_file = output_override;

        Model model = ReadMdpa(settings.mdpa_file);
        ApplyJsonProcesses(model, settings.project);

        if (!load_overrides.empty()) {
            for (auto& node : model.nodes) {
                node.loads = {};
            }
        }
        for (const auto& [node_id, load] : load_overrides) {
            Node& node = model.nodes.at(model.node_index_by_id.at(node_id));
            for (std::size_t d = 0; d < DofsPerNode; ++d) node.loads[d] += load[d];
        }
        for (const int node_id : fixed_node_overrides) {
            Node& node = model.nodes.at(model.node_index_by_id.at(node_id));
            node.fixed.fill(true);
        }
        for (const auto& [node_id, dof, value] : prescribed_overrides) {
            Node& node = model.nodes.at(model.node_index_by_id.at(node_id));
            node.fixed[dof] = true;
            node.values[dof] = value;
        }

        SparseRows K;
        Assemble(model, ReadMaterial(settings.materials_file, model), K);
        const std::vector<double> U = SolveReducedBandedCholesky(K, model);
        WriteVtk(settings.output_file, model, U);

        std::cout << "Read " << model.nodes.size() << " nodes, " << model.elements.size()
                  << " elements. Wrote " << settings.output_file << '\n';
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
