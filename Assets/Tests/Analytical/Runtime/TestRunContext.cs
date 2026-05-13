using UnityEngine;

namespace ShrodingerFlow.AnalyticalTests
{
    /// <summary>
    /// Контекст одного запуска теста: compute-шейдеры, каталог вывода и колбэки прогресса для GUI.
    /// Создаётся <see cref="AnalyticalTestRunner"/> и передаётся в <see cref="AnalyticalTestBase.RunTest"/>.
    /// </summary>
    public sealed class TestRunContext
    {
        public ComputeShader KernelsShader;
        public ComputeShader FftShader;
        public ComputeShader LesShader;
        /// <summary>Папка для артефактов теста (CSV / JSON / PNG). Создаётся раннером заранее.</summary>
        public string OutputDir;
        /// <summary>Камера для скриншотов (может быть null — тогда скриншот через ScreenCapture).</summary>
        public Camera Camera;
        /// <summary>Идентификатор запуска.</summary>
        public string RunId;

        // ---- Live GUI feedback ----

        /// <summary>Вызывается тестом для обновления прогресс-бара: (currentStep, totalSteps, statusLine).</summary>
        public System.Action<int, int, string> OnProgress;

        /// <summary>
        /// Вызывается тестом для передачи preview-текстуры (например срез |ω|).
        /// Раннер отображает её в Game View. Текстура не уничтожается контекстом — тест управляет временем жизни.
        /// </summary>
        public System.Action<Texture2D> OnPreviewTexture;

        /// <summary>Сообщить о текущем прогрессе.</summary>
        public void ReportProgress(int step, int totalSteps, string statusLine = null) =>
            OnProgress?.Invoke(step, totalSteps, statusLine);

        /// <summary>Передать preview-текстуру для отображения в GUI.</summary>
        public void SetPreview(Texture2D tex) =>
            OnPreviewTexture?.Invoke(tex);
    }
}
