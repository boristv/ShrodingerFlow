using System.Collections;
using UnityEngine;

namespace ShrodingerFlow.AnalyticalTests
{
    /// <summary>
    /// Базовый класс для одного аналитического теста: вешается рядом с <see cref="AnalyticalTestRunner"/>,
    /// раннер собирает все включённые компоненты и поочерёдно запускает <see cref="RunTest"/>.
    /// </summary>
    public abstract class AnalyticalTestBase : MonoBehaviour
    {
        [SerializeField] private bool _enabledInBatch = true;

        public bool IsEnabledInBatch => _enabledInBatch && enabled;

        /// <summary>Имя теста: используется как имя подпапки в каталоге результатов (без пробелов).</summary>
        public abstract string TestName { get; }

        /// <summary>Опциональное человекочитаемое описание для лога.</summary>
        public virtual string Description => string.Empty;

        /// <summary>
        /// Запускает тест. Должен быть корутиной (yield return null между шагами — иначе Unity заморозит фрейм).
        /// Все артефакты пишутся в <see cref="TestRunContext.OutputDir"/>.
        /// </summary>
        public abstract IEnumerator RunTest(TestRunContext ctx);
    }
}
