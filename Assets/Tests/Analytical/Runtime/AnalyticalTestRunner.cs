using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UnityEngine;
#if UNITY_EDITOR
using UnityEditor;
#endif

namespace ShrodingerFlow.AnalyticalTests
{
    /// <summary>
    /// Главный диспетчер аналитических тестов: держит ссылки на compute-шейдеры,
    /// собирает <see cref="AnalyticalTestBase"/> с того же GameObject и запускает их по очереди.
    /// Каждому тесту выдаётся персональный каталог <c>AnalyticalTestResults/&lt;TestName&gt;/&lt;runId&gt;/</c>.
    /// </summary>
    [DisallowMultipleComponent]
    public sealed class AnalyticalTestRunner : MonoBehaviour
    {
        [Header("Compute Shaders (ISF / FFT / LES)")]
        [SerializeField] private ComputeShader _kernelsShader;
        [SerializeField] private ComputeShader _fftShader;
        [SerializeField] private ComputeShader _lesShader;

        [Header("Output")]
        [Tooltip("Папка относительно корня репозитория (рядом с Assets/). Создаётся автоматически.")]
        [SerializeField] private string _outputRoot = "AnalyticalTestResults";

        [Header("Auto-run")]
        [Tooltip("Запустить все включённые тесты сразу после Play.")]
        [SerializeField] private bool _autoRunOnPlay = true;
        [Tooltip("Выйти из Play-режима после завершения всех тестов (удобно для batch-запуска).")]
        [SerializeField] private bool _stopPlayOnComplete = true;
        [Tooltip("Для скриншотов: камера, с которой снимать. Если не задана — используется Camera.main.")]
        [SerializeField] private Camera _screenshotCamera;

        [Header("Визуальная обратная связь (Game View)")]
        [Tooltip("Показывать оверлей с прогрессом и метриками.")]
        [SerializeField] private bool _showGui = true;
        [Tooltip("Показывать preview-текстуру (срез |ω| или ρ), если тест её передаёт.")]
        [SerializeField] private bool _showPreviewTexture = true;
        [Tooltip("Ширина GUI-панели в пикселях.")]
        [SerializeField] private float _guiPanelWidth = 460f;

        // ---- состояние GUI ----
        private bool _running;
        private string _guiCurrentTest = "";
        private int _guiStep;
        private int _guiTotalSteps;
        private string _guiStatusLine = "";
        private int _guiTestIndex;
        private int _guiTestTotal;
        private double _guiTestStartTime;
        private readonly System.Collections.Generic.List<string> _guiLog =
            new System.Collections.Generic.List<string>();
        private Texture2D _guiPreview;

        public string ResolvedOutputRoot => Path.Combine(GetRepoRoot(), _outputRoot);

        private void Start()
        {
            if (_autoRunOnPlay)
                StartCoroutine(RunAll());
        }

        private void OnGUI()
        {
            if (!_showGui || !Application.isPlaying) return;

            float panelW = _guiPanelWidth;
            float x = 10f;
            float y = 10f;
            float lineH = 22f;
            float pad = 8f;

            // Высота: заголовок + прогресс-бар + статус + лог (до 6 строк) + превью
            int logLines = Mathf.Min(_guiLog.Count, 6);
            float previewH = (_showPreviewTexture && _guiPreview != null) ? panelW * 0.5f + pad : 0f;
            float panelH = lineH * (3 + logLines) + 28f + pad * 2 + previewH;

            // Фон
            var bgRect = new Rect(x - pad, y - pad, panelW + pad * 2, panelH);
            GUI.color = new Color(0f, 0f, 0f, 0.72f);
            GUI.DrawTexture(bgRect, Texture2D.whiteTexture);
            GUI.color = Color.white;

            var labelStyle = new GUIStyle(GUI.skin.label)
            {
                fontSize = 13,
                normal = { textColor = Color.white },
                wordWrap = false
            };
            var boldStyle = new GUIStyle(labelStyle) { fontStyle = FontStyle.Bold };
            var dimStyle = new GUIStyle(labelStyle)
            {
                normal = { textColor = new Color(0.7f, 0.85f, 1f) },
                fontSize = 12
            };

            if (!_running)
            {
                GUI.Label(new Rect(x, y, panelW, lineH), "Analytical Tests — ожидание...", boldStyle);
                return;
            }

            // Строка 1 — тест и индекс
            GUI.Label(new Rect(x, y, panelW, lineH),
                $"[{_guiTestIndex}/{_guiTestTotal}] {_guiCurrentTest}", boldStyle);
            y += lineH;

            // Строка 2 — время
            double elapsed = Time.realtimeSinceStartupAsDouble - _guiTestStartTime;
            GUI.Label(new Rect(x, y, panelW, lineH),
                $"Elapsed: {elapsed:F1} s", dimStyle);
            y += lineH;

            // Прогресс-бар
            if (_guiTotalSteps > 0)
            {
                float progress = Mathf.Clamp01((float)_guiStep / _guiTotalSteps);
                var barBgRect = new Rect(x, y, panelW, 18f);
                GUI.color = new Color(0.2f, 0.2f, 0.2f);
                GUI.DrawTexture(barBgRect, Texture2D.whiteTexture);
                GUI.color = new Color(0.2f, 0.75f, 0.35f);
                GUI.DrawTexture(new Rect(x, y, panelW * progress, 18f), Texture2D.whiteTexture);
                GUI.color = Color.white;
                GUI.Label(new Rect(x + 4f, y, panelW, 18f),
                    $"step {_guiStep} / {_guiTotalSteps}  ({progress * 100f:F0} %)", dimStyle);
            }
            y += 24f;

            // Строка статуса
            if (!string.IsNullOrEmpty(_guiStatusLine))
            {
                GUI.Label(new Rect(x, y, panelW, lineH), _guiStatusLine, dimStyle);
                y += lineH;
            }

            // Лог (последние строки)
            int startLog = Mathf.Max(0, _guiLog.Count - 6);
            for (int i = startLog; i < _guiLog.Count; i++)
            {
                GUI.Label(new Rect(x, y, panelW, lineH), _guiLog[i], dimStyle);
                y += lineH;
            }

            // Preview-текстура
            if (_showPreviewTexture && _guiPreview != null)
            {
                y += pad;
                float texW = panelW;
                float texH = texW * _guiPreview.height / Mathf.Max(1, _guiPreview.width);
                GUI.DrawTexture(new Rect(x, y, texW, texH), _guiPreview, ScaleMode.ScaleToFit);
            }
        }

        [ContextMenu("Run All Enabled Tests")]
        public void RunAllFromMenu()
        {
            if (!Application.isPlaying)
            {
                Debug.LogWarning("[AnalyticalTestRunner] Запускайте из Play-режима (Unity Editor).");
                return;
            }
            if (_running)
            {
                Debug.LogWarning("[AnalyticalTestRunner] Уже выполняется. Повторный запуск проигнорирован.");
                return;
            }
            StartCoroutine(RunAll());
        }

        public IEnumerator RunAll()
        {
            if (_running)
                yield break;

            _running = true;
            if (!ValidateShaders())
            {
                _running = false;
                yield break;
            }

            var tests = GetComponents<AnalyticalTestBase>()
                .Where(t => t.IsEnabledInBatch)
                .ToList();
            if (tests.Count == 0)
            {
                Debug.LogWarning("[AnalyticalTestRunner] Нет включённых тестов на этом объекте.");
                _running = false;
                yield break;
            }

            string runId = DateTime.Now.ToString("yyyy-MM-dd_HH-mm-ss");
            string root = ResolvedOutputRoot;
            Directory.CreateDirectory(root);
            Debug.Log($"[AnalyticalTestRunner] Запуск {tests.Count} теста(ов). Каталог: {root}");
            Debug.Log($"[AnalyticalTestRunner] RunId={runId}");

            var summary = new BatchSummary
            {
                runId = runId,
                startedAt = DateTime.Now.ToString("o"),
                outputRoot = root,
                tests = new List<BatchTestEntry>()
            };

            _guiTestTotal = tests.Count;
            _guiPreview = null;

            for (int i = 0; i < tests.Count; i++)
            {
                var test = tests[i];
                string testDir = Path.Combine(root, test.TestName, runId);
                Directory.CreateDirectory(testDir);
                Debug.Log($"[AnalyticalTestRunner] [{i + 1}/{tests.Count}] {test.TestName} → {testDir}");

                _guiCurrentTest = test.TestName;
                _guiTestIndex = i + 1;
                _guiStep = 0;
                _guiTotalSteps = 0;
                _guiStatusLine = test.Description;
                _guiTestStartTime = Time.realtimeSinceStartupAsDouble;
                _guiPreview = null;

                var ctx = new TestRunContext
                {
                    KernelsShader = _kernelsShader,
                    FftShader = _fftShader,
                    LesShader = _lesShader,
                    OutputDir = testDir,
                    Camera = _screenshotCamera != null ? _screenshotCamera : Camera.main,
                    RunId = runId,
                    OnProgress = (step, total, status) =>
                    {
                        _guiStep = step;
                        _guiTotalSteps = total;
                        if (!string.IsNullOrEmpty(status))
                            _guiStatusLine = status;
                    },
                    OnPreviewTexture = tex => { _guiPreview = tex; },
                };

                double startTime = Time.realtimeSinceStartupAsDouble;
                Exception err = null;
                IEnumerator inner = null;
                try
                {
                    inner = test.RunTest(ctx);
                }
                catch (Exception ex)
                {
                    err = ex;
                }
                if (err == null && inner != null)
                {
                    yield return RunChild(inner, e => err = e);
                }
                double elapsed = Time.realtimeSinceStartupAsDouble - startTime;

                summary.tests.Add(new BatchTestEntry
                {
                    name = test.TestName,
                    outputDir = testDir,
                    description = test.Description,
                    elapsedSeconds = (float)elapsed,
                    error = err != null ? err.Message : null,
                });

                string resultLine = err != null
                    ? $"✗ {test.TestName}: {err.Message}"
                    : $"✓ {test.TestName} ({elapsed:F1} s)";
                _guiLog.Add(resultLine);

                if (err != null)
                    Debug.LogError($"[AnalyticalTestRunner] {test.TestName} упал: {err}");
                else
                    Debug.Log($"[AnalyticalTestRunner] {test.TestName} OK ({elapsed:F2} c)");
            }

            summary.finishedAt = DateTime.Now.ToString("o");
            string summaryPath = Path.Combine(root, $"session_{runId}.json");
            File.WriteAllText(summaryPath, JsonUtility.ToJson(summary, true));
            Debug.Log($"[AnalyticalTestRunner] Готово. Сводка: {summaryPath}");

            _running = false;

            if (_stopPlayOnComplete)
            {
#if UNITY_EDITOR
                EditorApplication.isPlaying = false;
#else
                Application.Quit();
#endif
            }
        }

        private static IEnumerator RunChild(IEnumerator inner, Action<Exception> onError)
        {
            while (true)
            {
                object current;
                try
                {
                    if (!inner.MoveNext())
                        yield break;
                    current = inner.Current;
                }
                catch (Exception ex)
                {
                    onError(ex);
                    yield break;
                }
                yield return current;
            }
        }

        private bool ValidateShaders()
        {
            bool ok = true;
            if (_kernelsShader == null) { Debug.LogError("[AnalyticalTestRunner] Не назначен KernelsShader (SFComputeKernels)."); ok = false; }
            if (_fftShader == null) { Debug.LogError("[AnalyticalTestRunner] Не назначен FftShader (SFComputeFFT)."); ok = false; }
            if (_lesShader == null) { Debug.LogError("[AnalyticalTestRunner] Не назначен LesShader (SFComputeLES)."); ok = false; }
            return ok;
        }

        private static string GetRepoRoot()
        {
            return Directory.GetParent(Application.dataPath).FullName;
        }

        [Serializable]
        private struct BatchTestEntry
        {
            public string name;
            public string outputDir;
            public string description;
            public float elapsedSeconds;
            public string error;
        }

        [Serializable]
        private struct BatchSummary
        {
            public string runId;
            public string startedAt;
            public string finishedAt;
            public string outputRoot;
            public List<BatchTestEntry> tests;
        }
    }
}
