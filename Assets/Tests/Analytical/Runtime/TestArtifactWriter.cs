using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;
using UnityEngine;

namespace ShrodingerFlow.AnalyticalTests
{
    /// <summary>
    /// Поток-помощник, который пишет CSV / JSON / PNG / лог для одного теста.
    /// CSV открывается лениво при первом <see cref="CsvLog.Row"/>, чтобы пустые файлы не плодились.
    /// </summary>
    public sealed class TestArtifactWriter : IDisposable
    {
        private readonly string _dir;
        private readonly string _logPath;
        private readonly List<CsvLog> _csvs = new List<CsvLog>();
        private StreamWriter _logWriter;

        public string Dir => _dir;

        public TestArtifactWriter(string outputDir)
        {
            _dir = outputDir;
            Directory.CreateDirectory(_dir);
            _logPath = Path.Combine(_dir, "run.log");
        }

        /// <summary>Создаёт или переоткрывает CSV-файл. Заголовок пишется только при первом Row().</summary>
        public CsvLog CreateCsv(string fileName, params string[] header)
        {
            string path = Path.Combine(_dir, fileName);
            var csv = new CsvLog(path, header);
            _csvs.Add(csv);
            return csv;
        }

        /// <summary>JSON через UnityEngine.JsonUtility. Объект должен быть [Serializable].</summary>
        public void Json(string fileName, object payload, bool pretty = true)
        {
            string path = Path.Combine(_dir, fileName);
            string json = JsonUtility.ToJson(payload, pretty);
            File.WriteAllText(path, json);
        }

        /// <summary>Сохраняет Texture2D в PNG (для срезов |ω|, плотности и т.п.).</summary>
        public string Png(string fileName, Texture2D tex)
        {
            string path = Path.Combine(_dir, fileName);
            byte[] bytes = tex.EncodeToPNG();
            File.WriteAllBytes(path, bytes);
            return path;
        }

        /// <summary>Скриншот текущего фрейма (через ScreenCapture, асинхронно — файл появится после следующего рендера).</summary>
        public string Screenshot(string fileName)
        {
            string path = Path.Combine(_dir, fileName);
            ScreenCapture.CaptureScreenshot(path);
            return path;
        }

        /// <summary>Скриншот с конкретной камеры — синхронный (Render-to-Texture).</summary>
        public string CameraScreenshot(string fileName, Camera camera, int width, int height)
        {
            if (camera == null)
                return Screenshot(fileName);

            var prevTarget = camera.targetTexture;
            var rt = RenderTexture.GetTemporary(width, height, 24, RenderTextureFormat.ARGB32);
            camera.targetTexture = rt;

            var prevActive = RenderTexture.active;
            RenderTexture.active = rt;
            camera.Render();

            var tex = new Texture2D(width, height, TextureFormat.RGB24, false);
            tex.ReadPixels(new Rect(0, 0, width, height), 0, 0);
            tex.Apply(false);

            RenderTexture.active = prevActive;
            camera.targetTexture = prevTarget;
            RenderTexture.ReleaseTemporary(rt);

            string path = Path.Combine(_dir, fileName);
            File.WriteAllBytes(path, tex.EncodeToPNG());
            UnityEngine.Object.Destroy(tex);
            return path;
        }

        public void Log(string message)
        {
            if (_logWriter == null)
                _logWriter = new StreamWriter(_logPath, append: false) { AutoFlush = true };
            _logWriter.WriteLine($"[{DateTime.Now:HH:mm:ss}] {message}");
            Debug.Log($"[AnalyticalTest] {message}");
        }

        public void Dispose()
        {
            foreach (var csv in _csvs)
                csv.Dispose();
            _csvs.Clear();
            _logWriter?.Dispose();
            _logWriter = null;
        }
    }

    /// <summary>CSV-логгер для метрик: одна строка на семпл.</summary>
    public sealed class CsvLog : IDisposable
    {
        private readonly string _path;
        private readonly string[] _header;
        private StreamWriter _writer;
        private int _expectedColumns;

        public string Path => _path;

        public CsvLog(string path, string[] header)
        {
            _path = path;
            _header = header;
            _expectedColumns = header.Length;
        }

        public void Row(params object[] values)
        {
            if (_writer == null)
            {
                _writer = new StreamWriter(_path, append: false) { AutoFlush = true };
                _writer.WriteLine(string.Join(",", _header));
            }
            if (values.Length != _expectedColumns)
            {
                Debug.LogWarning($"[CsvLog] {System.IO.Path.GetFileName(_path)}: ожидалось {_expectedColumns} колонок, получено {values.Length}.");
            }
            var sb = new StringBuilder(values.Length * 12);
            for (int i = 0; i < values.Length; i++)
            {
                if (i > 0) sb.Append(',');
                AppendCell(sb, values[i]);
            }
            _writer.WriteLine(sb.ToString());
        }

        private static void AppendCell(StringBuilder sb, object value)
        {
            switch (value)
            {
                case null: break;
                case float f:
                    sb.Append(f.ToString("R", CultureInfo.InvariantCulture));
                    break;
                case double d:
                    sb.Append(d.ToString("R", CultureInfo.InvariantCulture));
                    break;
                case int i:
                    sb.Append(i.ToString(CultureInfo.InvariantCulture));
                    break;
                case bool b:
                    sb.Append(b ? "1" : "0");
                    break;
                case string s:
                    if (s.IndexOfAny(new[] { ',', '"', '\n' }) >= 0)
                        sb.Append('"').Append(s.Replace("\"", "\"\"")).Append('"');
                    else
                        sb.Append(s);
                    break;
                default:
                    sb.Append(Convert.ToString(value, CultureInfo.InvariantCulture));
                    break;
            }
        }

        public void Dispose()
        {
            _writer?.Dispose();
            _writer = null;
        }
    }
}
